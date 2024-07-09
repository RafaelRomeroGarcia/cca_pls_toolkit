function plot_weight_stem(res, T, wfname, mod)
% plot_weight_stem
%
% # Syntax
%   plot_weight_stem(res, T, wfname, mod)
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
%   pass input weight and weight index as table
%   remove weights not in subset based on res.behav.weight.display
%   remove true weight normalization

% Open figure
if ~isempty(res.gen.figure.Position)
    figure('Position', res.gen.figure.Position);
else
    figure;
end

idout = false(numel(T.weight), 1);

% Load true weight if exists
if exist(res.simul.weight.file.(mod), 'file')    
    % Load true weight
    data = load(res.simul.weight.file.(mod));
    field = fieldnames(data);
    weight_true = reshape(data.(field{1}), [], 1);

    % Sort weight
    [weight_true, T.index] = sort(weight_true, 'descend');

    if numel(weight_true) > 100
        idout = randperm(numel(weight_true))' < numel(weight_true)-99;
        weight_true(idout) = []; 
    end
end

% Reorder table
T = T(T.index,:);

% Keep only subset of weights 
if strcmp(res.behav.weight.display, 'subset')
    T(~T.subset,:) = [];
end

% Normalize weight
if strcmp(res.simul.weight.norm, 'minmax')
    minmax = max(abs(T.weight));
    T.weight = T.weight / minmax;
elseif strcmp(res.simul.weight.norm, 'std')
    T.weight = T.weight / std(T.weight);
end

% Plot weight
stem(find(~idout), T.weight);
lg = {'weight'};

% Plot labels
xlabel(res.simul.xlabel);
ylabel(res.simul.ylabel);

if exist(res.simul.weight.file.(mod), 'file')    
    % Normalize weight
    if strcmp(res.simul.weight.norm, 'minmax')
        minmax = max(abs(weight_true));
        weight_true = weight_true / minmax;
    elseif strcmp(res.simul.weight.norm, 'std')
        weight_true = weight_true / std(weight_true);
    end

    % Plot true weight
    hold on;
    stem(find(~idout), weight_true);
    lg = [lg {'true weight'}];
end

% Update legend and axes
name_value = parse_struct(res.gen.legend);
legend(lg, name_value{:});
name_value = parse_struct(res.gen.axes);
parse_input(gca, name_value{:});

saveas(gcf, [wfname res.gen.figure.ext]);