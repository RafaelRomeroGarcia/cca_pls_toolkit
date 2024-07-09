function [param, S, bestid] = get_hyperparam(res, opt)
% get_hyperparam
%
% # Syntax
%   [param, S, bestid] = get_hyperparam(res, opt)
%
%_______________________________________________________________________
% Copyright (C) 2022 University College London
%
% Written by Agoston Mihalik (cca-pls-toolkit@cs.ucl.ac.uk)
% $Id$
%
% This file is part of CCA/PLS Toolkit.
%
% CCA/PLS Toolkit is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CCA/PLS Toolkit is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CCA/PLS Toolkit. If not, see <https://www.gnu.org/licenses/>.
%
% Modified by Agoston Mihalik (am3022@cantab.ac.uk)
%   fix grid search of matched hyperparameters
%   refactor code if opt is not default by introducing nreps

% Load cfg
if isfield(res.dir, 'frwork')
    cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');
else
    cfg = res;
end

% Number of hyperparameter levels
p = cfg.machine.param.name; % shorthand variable
num = zeros(1, numel(p));
for i=1:numel(p)
    num(i) = numel(cfg.machine.param.(p{i}));
end

% Compute design for hyperparameter combinations
if strcmp(cfg.machine.param.type, 'factorial')
    design = fullfact(num);
elseif strcmp(cfg.machine.param.type, 'matched')
    design = arrayfun(@(x) linspace(1, x, max(num)), num, 'un', 0);
    design = cat(1, design{:})';
end
        
if strcmp(opt, 'default')
    % Assign hyperparameters
    for i=1:size(design, 1)
        for j=1:size(design, 2)
            param(i).(p{j}) = cfg.machine.param.(p{j})(design(i,j));
        end
    end
        
else
    % Load default params
    param = get_hyperparam(res, 'default');
    nparams = numel(param);
    
    % Load compiled grid search file and format data
    if exist_file(cfg, fullfile(res.dir.(opt), 'allgrid.mat'))
        S = loadmat_struct(res, fullfile(res.dir.(opt), 'allgrid.mat'));
        for i=1:numel(cfg.machine.metric)
            nreps = size(S.(cfg.machine.metric{i}), 1) / nparams;
            S.(cfg.machine.metric{i}) = reshape(S.(cfg.machine.metric{i}), ...
                nreps, nparams, []);
        end
    else
        error('allgrid.mat not available');
    end

    % Calculate mean metric across subsamples
    fields = fieldnames(S);
    for i=1:numel(fields)
        S.(fields{i}) = real(nanmean(S.(fields{i}), 3));
    end
    
    % Find best hyperparameter                                          % deal with multiple values!!!
    bestid = zeros(1, nreps);
    switch cfg.machine.param.crit
        case 'correl+simwxy' % minimum distance
            for i=1:nreps
                [minval, bestid(i)] = min(calc_distance(S.correl(i,:), ...
                    nanmean([S.simwx(i,:); S.simwy(i,:)], 1)));
            end
            
        case 'correl' % maximum test correlation
            mx_crr = max(S.(cfg.machine.param.crit), [], 2);
            for i=1:nreps
                id = find(S.(cfg.machine.param.crit)(i,:) == mx_crr(i), 1); % sparsest solution if multiple values
                bestid(i) = id;
            end
    end
    param = param(bestid);
end
