function calc_bootid(res)
% calc_bootid
%
% Calculate sample indices for boostrapping analysis
%
% # Syntax:
%   calc_bootid(res)
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

% Load cfg
cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');

% Initialize file
file = fullfile(res.dir.frwork, 'boot', sprintf('bootmat_%d.mat', res.stat.nboot));

% Quit if bootstrap samples already generated
if res.stat.nboot == 0 || exist_file(res, file)
    return
end

% Set seed for reproducibility
rng(cfg.env.seed.boot);

% Load training and test indexes
[otrid, oteid] = loadmat(res, fullfile(res.dir.frwork, 'outmat.mat'), 'otrid', 'oteid');

% Calculate bootstrap samples
bootid = cell(1, res.frwork.split.nall);
for i=1:res.frwork.split.nall
    nsamples = sum(otrid(:,res.frwork.split.all(i)));
    bootid{i} = randi(nsamples, nsamples, res.stat.nboot);
end

% Save bootstrap samples
if ~isdir(fullfile(res.dir.frwork, 'boot'))
    mkdir(fullfile(res.dir.frwork, 'boot'))
end
savemat(res, file, 'bootid', bootid);