function res = save_results(res)
% save_results
%
% # Syntax
%   res = save_results(res)
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
%   update results table to use column headings based on cfg.env.save.tableHeading
%   update npcax/npcay column headings to be used independently in results table
%   use split/set as column heading in results table based on framework
%   fix varx/vary to be be used as table headings in results table
%   enable calculating number of PCs for summary table post-hoc

% Load cfg
cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');

% Load true metrics
S = loadmat_struct(res, fullfile(res.dir.res, 'model.mat'));

% Save res
if res.stat.sig
    % Significant results
    switch cfg.defl.crit
        case {'correl' 'covar'}
            save_res_pos(res, S.(cfg.defl.crit), {'max'});

        case 'pval+correl'
            save_res_pos(res, [res.stat.pval S.correl], {'min' 'max'});

        case 'correl+simwxy'
            distance = calc_distance(S.correl, nanmean(cat(2, S.simwx, S.simwy), 2));
            save_res_pos(res, distance, {'min'});
    end
else
    % No significant results
    save_res_neg(res);
end

% Write results table
output = {};
param = loadmat(res, fullfile(res.dir.res, 'param.mat'), 'param');
for i=1:numel(cfg.env.save.tableHeading)
    if contains(cfg.env.save.tableHeading{i}, {'nfeatx' 'nfeaty' 'npcax' 'npcay' ...
            'l2x' 'l2y' 'varx' 'vary'})
        m = cfg.env.save.tableHeading{i}(end); % shorthand variable
    end

    output = [output cfg.env.save.tableHeading(i)];

    switch cfg.env.save.tableHeading{i}
        case {'split' 'set'}
            output = [output {res.frwork.split.all}];

        case cfg.stat.crit
            output = [output {S.(cfg.stat.crit)}];

        case {'pval' 'correl'}
            output = [output {res.stat.pval}];

        case {'dist'}
            output = [output {distance}];

        case {'nfeatx' 'nfeaty'}
            output = [output {S.(['w' upper(m)])}];

        case {'l2x' 'l2y'}
            output = [output {cat(1, param.(['L2' m]))}];

        case {'varx' 'vary'}
            output = [output {cat(1, param.(['VAR' m]))}];

        case {'npcax' 'npcay'}
            if isfield(param, ['PCA' m])
                output = [output {cat(1, param.(['PCA' m]))}];
            else
                % Calculate number of PCs post-hoc
                npca = zeros(res.frwork.split.nall, 1);
                for j=1:res.frwork.split.nall
                    fname = getfname_dirload(cfg, 'svd', m, 'osplit', res.frwork.split.all(j));
                    if exist_file(cfg, fullfile(cfg.dir.load, 'svd', ['tr_' fname]))
                        trdata.(['L' m]) = loadmat(cfg, fullfile(cfg.dir.load, 'svd', ['tr_' fname]), ['L' upper(m)]);
                        featid = get_featid(trdata, param(j), m);
                        npca(j) = sum(featid);
                    end
                end
                output = [output {npca}];
            end

    end
end
write_results(res, 'results_table', output{:});


% --------------------------- Private functions ---------------------------

function write_results(res, fname, varargin)

% Create table
T = table();
for i=1:numel(varargin)
    if ~mod(i, 2)
        switch varargin{i-1}
            case {'split' 'set' 'npcax' 'npcay'}
                T.(varargin{i-1}) = varargin{i};

            case {'nfeatx' 'nfeaty'}
                T.(varargin{i-1}) = sum(varargin{i}~=0, 2);

            case {'dist' 'l2x' 'l2y' 'pval' 'correl' 'covar' 'varx' 'vary'}
                T.(varargin{i-1}) = arrayfun(@(x) sprintf('%.4f', x), ...
                    varargin{i}, 'un', 0);
        end
    end
end

% Write results
writetable(T, fullfile(res.dir.res, [fname '.txt']), 'Delimiter', '\t');


function res = save_res_pos(res, varargin)

% Assign input
[metric, fun] = varargin{1:2};
num = size(metric, 2);
if num ~= numel(fun)
    error('number of elements in metric and fun should match')
end

% Best split using criterions in specific order
bestid = (1:res.frwork.split.nall)';
for i=1:num
    fh = str2func(fun{i});
    met = cat(1, metric(bestid,i));
    [M, I] = fh(met);
    bestid = bestid(met == M); % we want all solutions
end
if numel(bestid) ~= 1
    warning('Multiple best splits, first index chosen');
end
res.frwork.split.best = res.frwork.split.all(bestid);

savemat(res, fullfile(res.dir.res, 'res.mat'), 'res', res);

% Display message based on verbosity level
switch res.env.verbose
    case {1 2}
        fprintf('\nSignificant results found!\n\n');
    case 3
        fprintf('Significant results found!\n\n');
end


function save_res_neg(res)

savemat(res, fullfile(res.dir.res, 'res.mat'), 'res', res);

% Display message based on verbosity level
switch res.env.verbose
    case {1 2}
        fprintf('\nNo significant results found at this level!\n\n');
    case 3
        fprintf('No significant results found at this level!\n\n');
end