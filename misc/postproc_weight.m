function T = postproc_weight(res, transform, T, varargin)
% postproc_weight
%
% # Syntax
%   T = postproc_weight(res, transform, T, varargin)
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
%   refactor code by introducing the following transforms: flip, positive, negative, top, minmax, sort, clean, strength, loading
%   add transform for calculating significant subset of weights based on bootstrapping analysis
%   add transform for calculating errorbars to weights based on bootstrapping analysis

% Load cfg
cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');

% Calculate loading
switch transform
    case 'flip'
        %----- Flip all weights to enable model comparison

        % Flip weights
        vars = T.Properties.VariableNames(ismember(T.Properties.VariableNames, {'weight' 'boot_weight'}));
        for i=1:numel(vars)
            T.(vars{i}) = -T.(vars{i});
        end

    case 'positive'
        %----- Filter weights by keeping positive values

        % Select subset
        T.subset = T.weight > 0;

    case 'negative'
        %----- Filter weights by keeping negative values

        % Select subset
        T.subset = T.weight < 0;

    case 'significant'
        %----- Filter weights by keeping only significant ones
        
        % Select subset
        if strcmp(res.gen.weight.stat.ci, 'sd')
            % Calculate standard deviation
            sd = nanstd(T.boot_weight, [], 2);

            % Parametric statistic based on scale of std
            T.subset = abs(T.weight ./ sd) > res.gen.weight.stat.sd;

        elseif strcmp(res.gen.weight.stat.ci, 'pi')
            % Calculate percentile interval
            pi = prctile(T.boot_weight, res.gen.weight.stat.pi, 2);

            % Non-parametric statistic based on width of percentile interval
            T.subset = min(pi .* sign(T.weight), [], 2) > 0;
        end

    case 'top'
        %----- Filter weights by keeping top ones based on absolute value
        
        % Top weights based on absolute value
        if isinf(res.gen.weight.subset.num)
            error('res.gen.weight.subset.num should be set for selecting top weights')
        end
        
        % Sort values
        [~, I] = sort(abs(T.(res.gen.weight.subset.order)), 'descend');
        
        % Select subset
        T.subset(I(res.gen.weight.subset.num+1:end)) = false;


    case 'minmax'
        %----- Filter weights by keeping top positive and negative ones

        % Top minimum and maximum weights
        if isinf(res.gen.weight.subset.num)
            error('res.gen.weight.subset.num should be set for selecting minmax weights')
        end

        % Sort values
        [~, I] = sort(T.(res.gen.weight.subset.order), 'descend');
        numnonneg = [sum(T.(res.gen.weight.subset.order)>0) sum(T.(res.gen.weight.subset.order)<0)];
        
        % Select subset
        if any(numnonneg >= res.gen.weight.subset.num)
            numnonneg(numnonneg>=res.gen.weight.subset.num) = res.gen.weight.subset.num;
            T.subset(I(numnonneg(1)+1:end-numnonneg(2))) = false;
        end

    case 'sort'
        %----- Sort weights for ordering labels on figure
        
        % Assign inputs
        modtype = varargin{:};
        
        % Sort values
        [~, T.index] = sort(T.(res.(modtype).label.order), 'descend');     
    
    case 'error'
        %----- Calculate errorbars for figure

        % Assign inputs
        modtype = varargin{:};

        if strcmp(res.(modtype).weight.errorbar.ci, 'sd')
            % Calculate scaled standard deviation 
            err = nanstd(T.boot_weight, [], 2) * res.(modtype).weight.errorbar.sd;

            % One or two-sided errorbar
            if strcmp(res.behav.weight.errorbar.side, 'one')
                T.errorbar(T.weight < 0, 1) =  err(T.weight < 0);
                T.errorbar(T.weight > 0, 2) =  err(T.weight > 0);
            elseif strcmp(res.behav.weight.errorbar.side, 'two')
                T.errorbar = repmat(err, 1, 2);
            end

        elseif strcmp(res.(modtype).weight.errorbar.ci, 'pi')
            % Calculate percentile
            pi = prctile(T.boot_weight, res.(modtype).weight.errorbar.pi, 2);
            
            % Two-sided errorbar
            if strcmp(res.behav.weight.errorbar.side, 'two')
                T.errorbar = [T.weight - pi(:,1) pi(:,2) - T.weight];
            else
                error('Percentile interval should be used with two-sided errorbars.')
            end
        end

    case 'clean'
        %----- Clean variables if they are a bit messy

        % Read label file
        labelfname = select_file(res, fullfile(res.dir.project, 'data'), ...
            'Select delimited label file for behaviour...', 'any', res.behav.file.label);
        T = cat(2, T, readlabel(labelfname, 'behavioural', T.Properties.VariableNames));

        % Flip specific weights for reversed-scored questionnaires
        if ismember('Flip', T.Properties.VariableNames)
            T.weight(T.Flip==1,1) = -T.weight(T.Flip==1,1);
        end

        % Set weight of secondary/highly redundant variables to 0 (e.g., see HCP)
        if ismember('Label_proc', T.Properties.VariableNames)
            [C, ia, ic] = unique(T.Label_proc);
            for i=1:numel(ia)
                iditem = find(ic==ia(i));
                if numel(iditem) > 1
                    w = abs(T.weight(iditem,1));
                    [M, I] = min(w);
                    T.subset(iditem(I),1) = false;
                end
            end
        end

    case 'strength'
        %----- Calculate strength by modifying weight by population mean data
        
        % Assign inputs
        mod = varargin{:};

        % Load original data
        data = load(res.data.(mod).fname);

        % Calculate strength
        vars = T.Properties.VariableNames(ismember(T.Properties.VariableNames, {'weight' 'boot_weight'}));
        for i=1:numel(vars)
            T.(vars{i}) = T.(vars{i}) .* reshape(sign(nanmean(data.(mod))), [], 1);
        end

    case 'loading'
        %----- Calculate loading by correlating input data and projection
        
        % Assign inputs
        [mod, split] = varargin{:};

        vars = T.Properties.VariableNames(ismember(T.Properties.VariableNames, {'weight' 'boot_weight'}));
        for i=1:numel(vars)
            % Assign weight
            weight = T.(vars{i});

            % Project data in feature space
            if ismember(cfg.machine.name, {'cca' 'rcca'})
                % Load parameters
                param = loadmat(res, fullfile(res.dir.res, 'param.mat'), 'param');
                param = param(split);

                % Load data in principal component basis (see transform_data.m)
                [trdata, trid, tedata, teid] = load_data(res, {['R' mod]}, ...
                    'osplit', split, param);
                data = concat_data(trdata, tedata, {['R' mod]}, trid, teid);

                % Define feature index
                featid = get_featid(trdata, param, mod);

                % RCCA inverse feature transformation
                weight = transform_data(cfg, 'svd-inverse-transform', weight, mod, ...
                    'V', trdata.(['V' mod])(:,featid));

                % Project data in feature space
                P = calc_proj(data.(['R' mod])(:,featid), weight);
            end

            % Data in input space
            [trdata, trid, tedata, teid] = load_data(res, {mod}, 'osplit', split);
            data = concat_data(trdata, tedata, {mod}, trid, teid);

            if ismember(cfg.machine.name, {'pls' 'spls'})
                % Project data in input space
                P = calc_proj(data.(mod), weight);
            end

            % Calculate correlation between input data and projection
            T.(vars{i}) = corr(P, data.(mod), 'rows', 'pairwise')';
        end
end