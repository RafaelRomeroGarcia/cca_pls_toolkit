function plot_weight_brain_node(res, T, wfname)
% plot_weight_brain_node
%
% # Syntax
%   plot_weight_brain_node(res, T, wfname)
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
%   pass input weight as table
%   add readlabel file to read connectivity label file

% Load label file
labelfname = select_file(res, fullfile(res.dir.project, 'data'), ...
    'Select delimited label file for regions...', 'any', res.roi.file.label);
T = cat(2, T, readlabel(labelfname, 'ROI', T.Properties.VariableNames, ...
    {'X' 'Y' 'Z' 'Label'}));
if ~ismember('Color', T.Properties.VariableNames)
    T.Color = ones(size(T, 1), 1);
end

% Weights
T.Size = T.weight;

% Remove regions we don't want to visualize
if ~isempty(res.roi.out)
    T(ismember(T.Index, res.roi.out),:) = [];
end

% Save results to text file
T = T(T.Size~=0,:);
[B, I] = sort(abs(T.Size), 'descend'); % sort nodes for easier readibility
T = T(I,:);
writetable(table(T.Label, T.Size, T.Color, 'VariableNames', {'Node' 'Weight' 'Color'}), [wfname '.csv']);

% Prepare BrainNet files (from scratch to avoid bad files)
T.Size = abs(T.Size / max(abs(T.Size)) * 2 * pi); % rescale to improve visualization
fname = init_brainnet(res, 'labelfname', labelfname, 'T', T);

% Plot brainnet
if exist(fname.options, 'file')
    BrainNet_MapCfg(fname.surf, fname.node, [wfname res.gen.figure.ext], fname.options);
else
    BrainNet_MapCfg(fname.surf, fname.node, [wfname res.gen.figure.ext]);
end