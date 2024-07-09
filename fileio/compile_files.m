function compile_files(res, files)
% compile_files
%
% # Syntax
%   compile_files(res, files)
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
%   update for bootstrapping analysis

% Load cfg
cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');

% Get option
[pathstr, name, ext] = fileparts(files{1});
opt = regexp(name, '^[a-z]+', 'match');
opt = opt{:};

% Initialization
if strcmp(opt, 'boot')
    S = cell2struct(cell(2, numel(files)), {'wX' 'wY'});
else
    S = cell2struct(cell(numel(cfg.machine.metric), numel(files)), cfg.machine.metric);
end

% Load files
for i=1:numel(files)
    if ~exist_file(res, fullfile(pathstr, ['all' opt '.mat']))
        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('%d\n', i);
            otherwise
                % display nothing at the moment
        end

        try
            S(i) = loadmat_struct(res, files{i});
        catch
            % Display message based on verbosity level
            switch res.env.verbose
                case 1
                    fprintf('%s was not loaded\n', files{i});
                otherwise
                    % display nothing at the moment
            end
        end
    else
        break
    end
end

% Compile files if no missing data
fname = fieldnames(S);
if strcmp(opt, 'boot') && ~any(arrayfun(@(x) isempty(S(x).wX), 1:numel(files)))
    name_value = parse_struct(S, 2);
    savemat(res, fullfile(pathstr, ['all' opt '.mat']), name_value{:});
    
elseif ~any(arrayfun(@(x) isempty(S(x).(fname{1})), 1:numel(files)))
    if strcmp(opt, 'grid')
        name_value = parse_struct(S, 1);
    elseif strcmp(opt, 'perm')
        name_value = parse_struct(S, 2);
    end
    savemat(res, fullfile(pathstr, ['all' opt '.mat']), name_value{:});
end
