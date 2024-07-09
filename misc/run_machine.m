function S = run_machine(cfg, trdata, tedata, featid, param, boot)
% run_machine
%
% # Syntax
%   S = run_machine(cfg, trdata, tedata, featid, param, boot)
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
%   use transform_data.m for SVD left-side feature transformations
%   update for bootstrapping analysis

% Initialize inputs
if ~exist('boot', 'var')
    boot = 0;
end

%----- Model training

% Run machine
switch cfg.machine.name
    case {'pls' 'spls'}
        % PLS/SPLS solved by power method
        [S.wX, S.wY] = spls(cfg, trdata, param);
        
        if boot
            % Deflate data and get all other components due to potential 
            % sign flipping and axis rotation in bootstrapping
            % ncomp = min([cfg.data.X.nfeat cfg.data.Y.nfeat]);
            ncomp = min([rank(trdata.X) rank(trdata.Y)]);

            for i=1:ncomp-1
                switch cfg.defl.name
                    case 'pls-projection'
                        for m='XY'
                            w = S.(['w' m])(:,i); % shorthand variable
                            
                            % Deflation step
                            trdata = deflation('pls-projection', trdata, m, w);
                        end
                        
                    case 'pls-modeA'
                        for m='XY'
                            w = S.(['w' m])(:,i); % shorthand variable

                            % Loading based on training data
                            p = trdata.(m)' * (trdata.(m) * w) / ((trdata.(m) * w)' * (trdata.(m) * w));

                            % Deflation step
                            trdata = deflation('pls-modeA', trdata, m, w, p);
                        end

                    case 'pls-regression'
                        w = S.wX(:,i); % shorthand variable

                        % Loading based on input training data
                        p = trdata.X' * (trdata.X * w) / ((trdata.X * w)' * (trdata.X * w));

                        % Deflation step
                        trdata = deflation('pls-regression', trdata, 'X', w, p);
                end

                % PLS/SPLS solved by power method
                [S.wX(:,i+1), S.wY(:,i+1)] = spls(cfg, trdata, param);
            end
        end
            
    case {'cca' 'rcca'}
        if ~boot
            % First component
            ncomp = 1;
        else
            % All components due to potential sign flipping and axis
            % rotation in bootstrapping
            ncomp = min([sum(featid.x) sum(featid.y)]);
        end

        % CCA/RCCA/PLS/PCA-CCA solved by standard eigenvalue problem
        [S.wX, S.wY] = rcca(trdata, featid, param, ncomp);
end

if ~boot
    %---- Model diagnostics

    % Compute projections for training data
    if ismember(cfg.machine.name, {'pls' 'spls'})
        trdata.PX = calc_proj(trdata.X, S.wX);
        trdata.PY = calc_proj(trdata.Y, S.wY);
    else
        trdata.PY = calc_proj(trdata.RY(:,featid.y), S.wY);
        trdata.PX = calc_proj(trdata.RX(:,featid.x), S.wX);
    end

    % Compute training metrics
    if ismember('trcorrel', cfg.machine.metric)
        % Correlation
        S.trcorrel = corr(trdata.PX, trdata.PY);
    end
    if ismember('trcovar', cfg.machine.metric)
        % Covariance
        S.trcovar = cov2(trdata.PX, trdata.PY);
    end
    isexvar = contains(cfg.machine.metric, 'trex');
    if any(isexvar)
        % Explained variance
        if strcmp(cfg.defl.name, 'generalized')
            S = calc_exvar(cfg, trdata, [], S, cfg.machine.metric(isexvar), param, featid);
        else
            S = calc_exvar(cfg, trdata, [], S, cfg.machine.metric(isexvar));
        end
    end
    cfg.machine.metric(isexvar) = [];

    %---- Model evaluation

    if ismember(cfg.machine.name, {'pls' 'spls'})
        tedata.PX = calc_proj(tedata.X, S.wX);
        tedata.PY = calc_proj(tedata.Y, S.wY);
    else
        tedata.PX = calc_proj(tedata.RX(:,featid.x), S.wX);
        tedata.PY = calc_proj(tedata.RY(:,featid.y), S.wY);
    end

    % Compute test metrics
    if ismember('correl', cfg.machine.metric)
        % Correlation
        S.correl = corr(tedata.PX, tedata.PY); % , 'Type', 'Spearman'
    end
    if ismember('covar', cfg.machine.metric)
        % Covariance
        S.covar = cov2(tedata.PX, tedata.PY);
    end
    isexvar = contains(cfg.machine.metric, 'ex');
    if any(isexvar)
        % Explained variance
        if strcmp(cfg.defl.name, 'generalized')
            S = calc_exvar(cfg, trdata, tedata, S, cfg.machine.metric(isexvar), param, featid);
        else
            S = calc_exvar(cfg, trdata, tedata, S, cfg.machine.metric(isexvar));
        end
    end
end

%---- Auxiliary steps

% Calculate primal weights
if contains(cfg.machine.name, 'cca')
    S.wX = transform_data(cfg, 'svd-left-transform', S.wX, 'X', 'V', trdata.VX(:,featid.x));
    S.wY = transform_data(cfg, 'svd-left-transform', S.wY, 'X', 'V', trdata.VY(:,featid.y));
end

if boot
    % Select weight most similar to true weight based on Procrustes analysis
    [S.wX, S.wY] = align_weight(trdata, S, 'select', cfg.machine.alignw);
end

% Record unsuccessful convergence for SPLS
if ismember('unsuc', cfg.machine.metric)
    S.unsuc = isnan(S.correl);
end