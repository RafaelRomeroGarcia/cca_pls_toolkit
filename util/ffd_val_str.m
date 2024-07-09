function str = ffd_val_str(cfg, varargin)
% ffd_val_str
%
% # Syntax
%   str = ffd_val_str(varargin)
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
%   fix for grid of matched hyperparameters
%   fix tabulation

% Parse input into field-value pairs
S = parse_input([], varargin{:});

% Number of values in fields
num = structfun(@numel, S);

if strcmp(cfg.machine.param.type, 'factorial')
    % Full factorial design of all value combinations
    dFF = fullfact(num);

elseif strcmp(cfg.machine.param.type, 'matched')
    % Matched design of hyperparameters
    param = {'L1x' 'L1y' 'L2x' 'L2y' 'PCAx' 'PCAy' 'VARx' 'VARy'}';
    num_param = num(ismember(fieldnames(S), param));
    design = arrayfun(@(x) linspace(1, x, max(num_param)), num_param, 'un', 0);
    design = cat(1, design{:})';

    % Full factorial design of restricted value combinations
    dFF = fullfact([num(~ismember(fieldnames(S), param)) size(design, 1)]);
    dFF = [dFF(:,1:end-1) design(dFF(:,end),:)];
end

% Create combinations in string format
str = cell(1, size(dFF, 1));
fn = fieldnames(S);
for i=1:size(dFF, 1)
    for j=1:size(dFF, 2)
        if mod(S.(fn{j})(dFF(i,j)), 1) == 0 % try decimal format
            str{i}{j} = [fn{j} '_' sprintf('%d', S.(fn{j})(dFF(i,j)))];

        elseif strfind(fn{j}, 'L2') % specific short format for L2 regularization
            str{i}{j} = [fn{j} '_' sprintf('%.4f', log(1-S.(fn{j})(dFF(i,j))))];

        else % otherwise use short format
            str{i}{j} = [fn{j} '_' sprintf('%.4f', S.(fn{j})(dFF(i,j)))];
        end
    end
    str{i} = strjoin(str{i}, '_');
end