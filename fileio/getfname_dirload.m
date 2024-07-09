function fname = getfname_dirload(cfg, op, m, splittype, split, varargin)
% getfname_dirload
%
% Get file name of calculations saved in the cfg.dir.load folder such as preproc and svd
%
% # Syntax
%   fname = getfname_dirload(cfg, op, m, splittype, split, varargin)
%___________________________________________________________________
% Copyright (C) 2023 MQ: Transforming Mental Health
%
% Written by Agoston Mihalik (am3022@cantab.ac.uk)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

switch op
    case 'preproc'
        % Create folder for preprocessing params
        if ~isdir(fullfile(cfg.dir.load, 'preproc'))
            mkdir(fullfile(cfg.dir.load, 'preproc'))
        end

        % Initialize file name
        if ismember(splittype, {'osplit' 'otrid' 'oteid' 'iext'})
            fname = sprintf('preproc%s_split_%d.mat', m, split(1));
        elseif ismember(splittype, {'isplit' 'itrid' 'iteid'})
            fname = sprintf('preproc%s_split_%d_subsample_%d.mat', m, split(1), split(2));
        elseif strcmp(splittype, 'oext')
            fname = sprintf('preproc%s_oext.mat', m);
        end


    case 'svd'
        % Create folder for svd results
        if ~isdir(fullfile(cfg.dir.load, 'svd'))
            mkdir(fullfile(cfg.dir.load, 'svd'))
        end

        % Initialize file name
        if ismember(splittype, {'osplit' 'otrid' 'oteid' 'iext'})
            fname = sprintf('svd%s_split_%d.mat', m, split(1));
        elseif ismember(splittype, {'isplit' 'itrid' 'iteid'})
            fname = sprintf('svd%s_split_%d_subsample_%d.mat', m, split(1), split(2));
        end
end
