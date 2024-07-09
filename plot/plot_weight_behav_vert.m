function plot_weight_behav_vert(res, T, wfname)
% plot_weight_behav_vert
%
% # Syntax
%   plot_weight_behav_vert(res, weight, iweight, wfname)
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
%   fix numeric category level in label file
%   pass input weight and weight index as table
%   add readlabel file to read behavioural label file
%   add color to weight table to enable explicit weight colors
%   remove or zero out weights not in subset based on res.behav.weight.display
%   add errorbars to weights based on bootstrapping analysis
%   refactor plotting bars to enable new functionalities above
%   add lines to delineate categories in plot if set by res.behav.categ.delineate
%   add xtick labels explicitly to plot

% Load label file
labelfname = select_file(res, fullfile(res.dir.project, 'data'), ...
    'Select delimited label file for behaviour...', 'any', res.behav.file.label);
T = cat(2, T, readlabel(labelfname, 'behavioural', T.Properties.VariableNames, {'Label'}));
iscategory = ismember('Category', T.Properties.VariableNames);
if ~iscategory
    T.Category = ones(size(T, 1), 1);
end
if isnumeric(T.Category)
    if res.behav.categ.delineate
        categ_line = find(diff(T.Category));
    end
    ndigits = arrayfun(@(x) numel(num2str(x)), T.Category);
    T.Category = sprintfc(['%0' num2str(max(ndigits)) 'd'], T.Category);
end

% Open figure
if ~isempty(res.gen.figure.Position)
    figure('Position', res.gen.figure.Position);
else
    figure;
end

% Normalize weight
if strcmp(res.behav.weight.norm, 'minmax')
    minmax = max(abs(T.weight));
    T.weight = T.weight / minmax;
elseif strcmp(res.behav.weight.norm, 'std')
    T.weight = T.weight / std(T.weight);
elseif strcmp(res.behav.weight.norm, 'zscore')
    T.weight = zscore(T.weight);
end

% Reorder table
T = T(T.index,:);

% Remove/zero out weights not in subset 
if strcmp(res.behav.weight.display, 'subset-all')
    T.weight(~T.subset) = 0;
elseif strcmp(res.behav.weight.display, 'subset-only')
    T(~T.subset,:) = [];
end
nweights = size(T, 1);

% Set colors for the categories
[categ, ia, ic] = unique(T.Category);
cmap = colormap('jet');
if size(cmap, 1) < numel(categ)
    error('Too many groups, not enough colors to plot them.')
end
cmap = cmap(round(linspace(1,size(cmap, 1),numel(categ))),:);
T.Color = cmap(ic,:);

% Plot weights
hold on;
x = 1:nweights;
h = zeros(1, nweights);
for i=1:nweights
    h(i) = bar(x(i), T.weight(i), 'FaceColor', T.Color(i,:));
end

% Plot error bar
if sum(sum(isnan(T.errorbar))) ~= size(T, 1) * 2
    id = T.weight ~= 0; % remove 0 weights
    er = errorbar(x(id), T.weight(id), T.errorbar(id, 1), T.errorbar(id, 2), 'vertical');
    set(er, 'Color', res.behav.weight.errorbar.Color, 'LineStyle', 'none', ...
        'LineWidth', res.behav.weight.errorbar.LineWidth);
end

% Add line to delineate categories
if res.behav.categ.delineate && i < numel(categ)
    plot(repmat(categ_line(i), 2, 1), [min(T.Weight); max(T.Weight)], ...
        'k--', 'LineWidth', 1.5);
end
hold off;

% Plot labels
ylabel(res.behav.ylabel);
xlabel(res.behav.xlabel);

% Update legend and axes
if numel(categ) > 1
    name_value = parse_struct(res.gen.legend);
    legend(h(ia), categ, name_value{:});
end
name_value = parse_struct(res.gen.axes);
set(gca, 'xTick', 1:numel(T.Label), 'xTickLabel', T.Label, name_value{:});
        
% Save figure
saveas(gcf, [wfname res.gen.figure.ext]);

% Save weights to csv
T.Properties.VariableNames{ismember(T.Properties.VariableNames, 'weight')} = 'Weight';
if iscategory
    writetable(T(:,{'Category' 'Label' 'Weight'}), [wfname '.csv'], ...
        'QuoteStrings', true);
else
    writetable(T(:,{'Label' 'Weight'}), [wfname '.csv'], 'QuoteStrings', true);
end