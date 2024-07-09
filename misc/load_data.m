function varargout = load_data(res, mod, splittype, split, param, bootid)
% load_data
%
% # Syntax
%   varargout = load_data(res, mod, splittype, split, param, bootid)
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
%   use transform_data.m for SVD right-side feature transformation
%   update for bootstrapping analysis
%   use getfname_dirload to get/generate filename in load directory
%   refactor code a bit   

% Load cfg
cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');

% Load training and test indexes
switch splittype
    case 'otrid'
        otrid = loadmat(res, fullfile(res.dir.frwork, 'outmat.mat'), 'otrid');
        trid = otrid(:,split);
        
    case 'itrid'
        itrid = loadmat(res, fullfile(res.dir.frwork, 'inmat.mat'), 'itrid', 'iteid');
        trid = itrid(:,split(1),split(2));
        
    case 'oext'
        oextid = loadmat(res, fullfile(res.dir.frwork, 'extmat.mat'), 'oextid');
        trid = oextid(:,split);
        
    case {'osplit' 'oteid'}
        [otrid, oteid] = loadmat(res, fullfile(res.dir.frwork, 'outmat.mat'), 'otrid', 'oteid');
        trid = otrid(:,split);
        teid = oteid(:,split);
        
    case {'isplit' 'iteid'}
        [itrid, iteid] = loadmat(res, fullfile(res.dir.frwork, 'inmat.mat'), 'itrid', 'iteid');
        trid = itrid(:,split(1),split(2));
        teid = iteid(:,split(1),split(2));
        
    case 'iext'
        otrid = loadmat(res, fullfile(res.dir.frwork, 'outmat.mat'), 'otrid');
        oextid = loadmat(res, fullfile(res.dir.frwork, 'extmat.mat'), 'oextid');
        
        % Concatenate original training indexes and external indexes
        trid = cat(1, otrid(:,split(1)), false(size(oextid, 1), 1));
        teid = cat(1, false(size(otrid, 1), 1), oextid(:,split(2)));
end

% Initialize data
[data, trdata, tedata] = deal(struct());

% Add confounds to modality if needed
if cfg.data.conf
    mod = [{'C'} mod];
end

for i=1:numel(mod)
    % Shorthand for modality without processing
    m = erase(mod{i}, 'R');
    lm = lower(m);

    proc = 0;
    if ~strcmp(mod{i}, 'C') && exist('param', 'var') && contains(cfg.machine.name, 'cca')
        % Initialize folder/file for SVD results
        fname_svd = getfname_dirload(cfg, 'svd', lm, splittype, split);
        if ~exist('bootid', 'var') && exist_file(cfg, fullfile(cfg.dir.load, 'svd', ['tr_' fname_svd])) ...
                && exist_file(cfg, fullfile(cfg.dir.load, 'svd', ['te_' fname_svd]))
            proc = 1;
        end
    end
    if ~exist('param', 'var') && ~strncmp(mod{i}, 'R', 1)
       proc = 0; 
    end
    
    if proc
        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('Loading initial SVD results of %s...', m);
                tic;
            otherwise
                % display nothing at the moment
        end

        % Load processed data
        [trdata.(['R' m]), trdata.(['V' m]), trdata.(['L' m])] = loadmat(cfg, ...
            fullfile(cfg.dir.load, 'svd', ['tr_' fname_svd]), ['R' m], ['V' m], ['L' m]);
        if exist('teid', 'var')
            tedata.(['R' m]) = loadmat(cfg, fullfile(cfg.dir.load, 'svd', ...
                ['te_' fname_svd]), ['R' m]);
        end

        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('done!\n');
                toc
            otherwise
                % display nothing at the moment
        end
    
    else
        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('Loading data %s...', m);
                tic;
            otherwise
                % display nothing at the moment
        end
        
        
        % Load original and/or external data
        switch splittype
            case 'oext'
                data = load_data_mod(res, data, m);
                
            case 'iext'
                data = load_data_mod(cfg, data, m);
                data = load_data_mod(res, data, m);
                
            otherwise
                data = load_data_mod(cfg, data, m);
        end
        
        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('done!\n');
                toc
            otherwise
                % display nothing at the moment
        end

        % Impute missing data
        data = impute_mat(cfg, data, trid, m);
        
        % Add bias term to confounds if needed
        if strcmp(m, 'C')
            if isfield(data, 'C') && ~any(arrayfun(@(x) isequal(ones(size(data.C, 1), 1), data.C(:,x)), 1:size(data.C, 2)))
                data.C = [ones(size(data.C, 1), 1) data.C];
            end
        end
        
        % Training examples in input space
        trdata.(m) = data.(m)(trid,:);
        
        % Test examples in input space
        if exist('teid', 'var')
            tedata.(m) = data.(m)(teid,:);
        end
        
        % Nothing else to do with confounds
        if strcmp(m, 'C')
            continue
        end
        
        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('Preprocessing data %s...', m);
            otherwise
                % display nothing at the moment
        end
        
        % Boostrap training data if needed
        if exist('bootid', 'var')
            trdata.(m) = trdata.(m)(bootid,:);
        end

        % Initialize folder/file for preprocessing params
        fname_preproc = getfname_dirload(cfg, 'preproc', lm, splittype, split);
        
        % Preprocess training and test data
        if exist('bootid', 'var')
            trdata = preproc_data(res, trdata, m, ['tr_boot_' fname_preproc]);
            if exist('teid', 'var')
                tedata = preproc_data(res, tedata, m, ['te_boot_' fname_preproc]);
            end
        else
            trdata = preproc_data(res, trdata, m, ['tr_' fname_preproc]);
            if exist('teid', 'var')
                tedata = preproc_data(res, tedata, m, ['te_' fname_preproc]);
            end
        end

        % Display progress based on verbosity level
        switch res.env.verbose
            case 1
                fprintf('done!\n');
                toc
            otherwise
                % display nothing at the moment
        end
        
        % Transform features
        if ismember(cfg.machine.name, {'pls' 'spls'}) ...
                || (~exist('param', 'var') && ~strncmp(mod{i}, 'R', 1))
            % No need for feature transformation here
        
        else
            % Display progress based on verbosity level
            switch res.env.verbose
                case 1
                    fprintf('Compute initial SVD for RAM/time-efficiency...');
                    toc
                otherwise
                    % display nothing at the moment
            end

            % Transform features into principal component basis
            if exist('bootid', 'var')
                trdata = transform_data(cfg, 'svd-right-transform', trdata, m, ...
                    'fname', ['boot_' fname_svd]);
                if exist('teid', 'var')
                    tedata = transform_data(cfg, 'svd-right-transform', tedata, m, ...
                        'V', trdata.(['V' m]));
                end

            else
                trdata = transform_data(cfg, 'svd-right-transform', trdata, m, ...
                    'fname', ['tr_' fname_svd]);
                if exist('teid', 'var')
                    tedata = transform_data(cfg, 'svd-right-transform', tedata, m, ...
                        'fname', ['te_' fname_svd], 'V', trdata.(['V' m]));
                end
            end

            % Display progress based on verbosity level
            switch res.env.verbose
                case 1
                    fprintf('done!\n');
                    toc
                otherwise
                    % display nothing at the moment
            end
        end
    end
    
    % Define feature index
    if exist('param', 'var')
        if ismember(cfg.machine.name, {'pls' 'spls'})
            featid.(lm) = true(1, cfg.data.(m).nfeat);
        else
            featid.(lm) = get_featid(trdata, param, m);
        end
    end
    
    % Deflate data
    if res.frwork.level > 1 ...
            && (isfield(trdata, ['R' m]) || ismember(cfg.machine.name, {'pls' 'spls'}))
        [trdata, tedata] = deflate_data(res, trdata, tedata, m, split(1));
    end
    
end

% if exist('param', 'var') && ismember(cfg.machine.name, {'pls' 'spls'})    
%     % Initialize folder/file for SVD results
%     if ~isdir(fullfile(cfg.dir.load, 'svd', sprintf('level%d', res.frwork.level)))
%         mkdir(fullfile(cfg.dir.load, 'svd', sprintf('level%d', res.frwork.level)))
%     end
%     if numel(split) == 1
%         fname = sprintf('tr_svdxy_split_%d.mat', split(1));
%     elseif numel(split) == 2
%         fname = sprintf('tr_svdxy_split_%d_subsample_%d.mat', split(1), split(2));
%     end
%     fname = fullfile(cfg.dir.load, 'svd', sprintf('level%d', res.frwork.level), fname);
%     
%     if ~exist_file(cfg, fname)
%         fprintf('Compute initial SVD for RAM/time-efficiency...');
%         
%         % Calculate cross-covariance matrix
%         trdata.XY = trdata.X' * trdata.Y;
%         
%         % Calculate SVD on cross-covariance matrix
%         trdata.VXY = fastsvd(trdata.XY, 1, 0, 1, 'V');
%         
%         % Save SVD results
%         name_value = parse_struct(trdata, 1, {'XY' 'VXY'});
%         savemat(res, fname, name_value{:});
%         
%     else
%         fprintf('Loading initial SVD results of XY...');
%         
%         % Load SVD results
%         [trdata.XY, trdata.VXY] = loadmat(res, fname, 'XY', 'VXY');
%     end
%     
%     fprintf('done!\n'); toc
% end

% Remove confounds from modality if needed
if isfield(cfg.data, 'C')
    mod(ismember(mod, 'C')) = [];
end

% Assign output in data+id format
varargout = {trdata trid};
switch splittype
    case {'oteid' 'iteid' 'osplit' 'isplit'}
        varargout = [varargout {tedata teid}];
        if exist('param', 'var')
            varargout{end+1} = featid;
        end
        
    case 'iext'
        varargout = [varargout {tedata teid(size(otrid, 1)+1:end)}];
end


% --------------------------- Private functions ---------------------------

function data = load_data_mod(res, data, mod)

% Load data modality
tmp = load(res.data.(mod).fname);

if ~isfield(data, mod)
    % Assign data
    data.(mod) = tmp.(mod);
else
    % Concatenate data - features should match, otherwise padded with 0!!
    [nsubj1, nfeat1] = size(data.(mod));
    [nsubj2, nfeat2] = size(tmp.(mod));
    if nfeat1 > nfeat2
        data.(mod) = cat(1, data.(mod), [tmp.(mod) zeros(nsubj2, nfeat1-nfeat2)]);
        warning('Testing features of %s padded with 0 to be able to concatenate data.', mod)
    elseif nfeat1 < nfeat2
        data.(mod) = cat(1, [data.(mod) zeros(nsubj1, nfeat2-nfeat1)], tmp.(mod));
        warning('Training features of %s padded with 0 to be able to concatenate data.', mod)
    else
        data.(mod) = cat(1, data.(mod), tmp.(mod));
    end
end
