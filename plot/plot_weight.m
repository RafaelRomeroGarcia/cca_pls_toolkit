function plot_weight(res, mod, modtype, split, func, varargin)
% plot_weight
%
% It plots the model weights in specific figures based on the modality of 
% the data. 
%
% # Syntax
%   plot_weight(res, mod, modtype, split, func, varargin)
%
% # Inputs
% res:: struct 
%   res structure containing information about results and plot specifications
% mod:: 'X', 'Y' 
%   modality of data to be used for plotting
% modtype:: 'behav', 'conn', 'vbm', 'roi', 'simul'
%   type of data
% split:: int
%   index of data split to be used
% func:: 'behav_horz', 'behav_vert', 'behav_text', 'brain_conn_node', 'brain_cortex', 'brain_edge', 'brain_module', 'brain_node', 'stem'
%   name of the specific plotting function (after `plot_weight_*` prefix) to
%   be called
% varargin:: name-value pairs
%   additional options can be passed via name-value pairs with dot notation
%   supported (e.g., 'behav.weight.numtop', 20)
%
% # Examples
%
% ## Modality independent
%   % Plot Y weights as stem plot
%   plot_weight(res, 'Y', 'simul', res.frwork.split.best, 'stem', ...
%   'gen.axes.YLim', [-0.2 1.2], 'simul.weight.norm', 'minmax', ...
%   'gen.axes.FontSize', 20, 'gen.legend.FontSize', 20);
%
% ![weight_stem](../figures/visualization_stem.png)
%
% ## Behaviour
%   % Plot behavioural weights as vertical bar plot
%   plot_weight(res, 'Y', 'behav', res.frwork.split.best, 'behav_vert', ...
%   'gen.axes.FontSize', 20, 'gen.legend.FontSize', 20, ...
%   'gen.axes.YLim', [-0.4 1.2], 'gen.weight.flip', 1, ...
%   'behav.weight.sorttype', 'sign', 'behav.weight.numtop', 20, ...
%   'behav.weight.norm', 'minmax');
%
% ![weight_behav](../figures/visualization_behav_vert.png)
%
% ## ROI-wise sMRI
%
%   % Plot ROI weights on a glass brain
%   plot_weight(res, 'X', 'roi', 1, 'brain_node', ...
%   'roi.weight.sorttype', 'sign', 'roi.weight.numtop', 20, ...
%   'roi.out', 9000 + [reshape([1:10:81; 2:10:82], [], 1); ...
%   reshape(100:10:170, [], 1)]);
%
% ![weight_roi](../figures/visualization_brain_node.png)
%
% ## fMRI connectivity edges
%   % Plot connectivity weights on glass brain
%   plot_weight(res, 'X', 'conn', res.frwork.split.best, 'brain_edge', ...
%   'conn.weight.sorttype', 'sign', 'conn.weight.numtop', 20);
%
% ![weight_conn_edge](../figures/visualization_brain_edge.png)
%
% ---
% See also: [plot_paropt](../plot_paropt), [plot_proj](../plot_proj)
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
%   remove loading cfg
%   move all weight postprocessing (e.g., calculating error bars, subseting and sorting weights) to postproc_weight.m
%   pass weight as table for plotting functions
%   add bootstrapping analysis
%   add calculation for significant subset of weights based on boostrapping analysis
%   add calculation for errorbars to weights based on bootstrapping analysis

% Parse input and add default settings
res = res_defaults(res, modtype, varargin{:});

% Add SPM if needed
if strcmp(res.gen.selectfile, 'interactive') || strcmp(modtype, 'vbm')
    set_path('spm');
end

%----- Get weight vectors

% Load weights
weight = loadmat(res, fullfile(res.dir.res, 'model.mat'), ['w' mod]);
T = table(weight(split,:)', 'VariableNames', {'weight'});

% Add boostrapped weights if exist and calculate Z statistic
if res.stat.nboot > 0
    % Load bootstrapped weights
    S_boot = loadmat_struct(res, fullfile(res.dir.boot, 'allboot.mat'));
    T.boot_weight = reshape(squeeze(S_boot.(['w' mod])(split,:,:))', [], res.stat.nboot);
    T.zstat = T.weight ./ nanstd(T.boot_weight, [], 2);
end

% Flip weights if needed
if res.gen.weight.flip
    T = postproc_weight(res, 'flip', T);
end

% Select subset of weights if needed
T.subset = true(size(T, 1), 1);
if ~strcmp(res.gen.weight.subset.type, 'none')
    T = postproc_weight(res, res.gen.weight.subset.type, T, modtype);
end

% Calculate loading if needed
if strcmp(res.gen.weight.type, 'correlation')    
    T = postproc_weight(res, 'loading', T, mod, split);
end

% Calculate strength if needed
if isfield(res.(modtype), 'weight') && isfield(res.(modtype).weight, 'type') ...
        && strcmp(res.(modtype).weight.type, 'strength')
    T = postproc_weight(res, 'strength', T, mod);
end
T.subset(isnan(T.weight)) = false; % not sure anymore why we need this but kept it to be safe

% Clean variables if needed
if isfield(res.(modtype), 'weight') && isfield(res.(modtype).weight, 'clean') ...
        && res.(modtype).weight.clean
    T = postproc_weight(res, 'clean', T);
end

% Sort weights if needed
T.index = (1:size(T, 1))';
if isfield(res.(modtype), 'label') && isfield(res.(modtype).label, 'order') ...
        && ~strcmp(res.(modtype).label.order, 'index')
    T = postproc_weight(res, 'sort', T, modtype);
end

% Calculate errorbars for weights if needed
T.errorbar = NaN(size(T, 1), 2);
if isfield(res.(modtype), 'weight') && isfield(res.(modtype).weight, 'errorbar') ...
        && ~strcmp(res.(modtype).weight.errorbar.ci, 'none') && res.stat.nboot > 0
    T = postproc_weight(res, 'error', T, modtype);
end

% Remove unnecessary variables from table
T(:,T.Properties.VariableNames(ismember(T.Properties.VariableNames, ...
    {'zstat' 'boot_weight'}))) = [];

%----- Visualize/summarize weights and write to disc

% Specify weight file name
if isfield(res.(modtype), 'weight') && isfield(res.(modtype).weight, 'type') ...
        && strcmp(res.(modtype).weight.type, 'strength')
    wfname = fullfile(res.dir.res, [res.gen.weight.type '_strength' mod]);
else
    wfname = fullfile(res.dir.res, [res.gen.weight.type mod]);
end
if strcmp(modtype, 'conn')
    suffix = strsplit(func, '_');
    wfname = [wfname '_' suffix{end}];
end
switch res.gen.weight.subset.type
    case {'positive' 'negative'}
        wfname = sprintf('%s_%s', wfname, res.gen.weight.subset.type(1:3));
    case 'significant'
        wfname = sprintf('%s_sig%s', wfname, res.gen.weight.stat.ci);
    case {'top' 'minmax'}
        wfname = sprintf('%s_%s%d', wfname, res.gen.weight.subset.type, res.gen.weight.subset.num);
end
if isfield(res.(modtype), 'module') && isfield(res.(modtype).module, 'type') ...
        && ~isempty(strfind(func, 'module'))
    wfname = [wfname '_' res.(modtype).module.type];
end
if isfield(res.(modtype), 'module') && isfield(res.(modtype).module, 'logtrans') ...
        && res.(modtype).module.logtrans
    wfname = [wfname '_log'];
end
wfname = sprintf('%s_split%d', wfname, split);

% Visualize/summarize weights
plottype = strsplit(func, '_');
plottype = plottype{1}; % first substring codes main plot type
func = str2func(['plot_weight_' func]);
if strcmp(plottype, 'behav')
    % Bar plots or text prints
    func(res, T, wfname);
    
elseif strcmp(plottype, 'stem')
    % Stem plot
    func(res, T, wfname, mod);
    
elseif strcmp(plottype, 'brain')
    % Glass brain plots (and module connectivity plot)
    func(res, T, wfname);
end