function [wX, wY] = align_weight(Sref, S, op, wm)      
% align_weight
%
% Align weight to a reference weight based on Procrustes analysis
%
% # Syntax:
%   [wX, wY] = align_weight(Sref, S, op, wm)
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

switch op
    case 'select'
        % Select bootstrapped weight most similar to true weight due to 
        % potential rotation and sign flipping
        sim = norm_features(Sref.(wm))' * norm_features(S.(wm));
        [~, iweight] = max(abs(sim));
        if sim(iweight) > 0
            wX = S.wX(:,iweight);
            wY = S.wY(:,iweight);
        else
            wX = -S.wX(:,iweight);
            wY = -S.wY(:,iweight);
        end

    case 'flip'
        % Flip weight to match reference weight
        if norm_features(Sref.(wm))' * norm_features(S.(wm)) < 0
            wX = -S.wX;
            wY = -S.wY;
        else
            wX = S.wX;
            wY = S.wY;
        end
end