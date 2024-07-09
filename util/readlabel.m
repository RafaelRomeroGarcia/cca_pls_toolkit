function T = readlabel(fname, labeltype, varsno, varsofint)
% readlabel
%
% # Syntax:
%   T = readlabel(fname, labeltype, varsno, varsofint)
%
%_______________________________________________________________________
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

T = readtable(fname);

% Check that variables in label file do not clash with existing variables
if ~isempty(intersect(varsno, T.Properties.VariableNames))
    error('The %s label file should not contain the following columns: %s', labeltype, strjoin(varsno, ', '));
end

% Check that label file contains all variables of interest
if nargin > 3
    if ~all(ismember(varsofint, T.Properties.VariableNames))
        error('The %s label file should contain the following columns: %s', labeltype, strjoin(varsofint, ', '));
    end
end