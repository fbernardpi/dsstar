%
% Simple demo code for the paper 
%
% DS*: Tighter Lifting-Free Convex Relaxations for Quadratic Matching Problems
% F. Bernard, C. Theobalt, M. Moeller
% IEEE Conference on Computer Vision and Pattern Recognition (CVPR). June, 2018
%
% We want to solve the quadratic assignment problem
% min_X X(:)'*W*X(:) s.t. X is a permutation matrix
%
%
% Author & Copyright (C) 2020: Florian Bernard (f.bernardpi[at]gmail[dot]com)
% 
% 
% LICENSE:
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

N = 15; % number of nodes that are matched
W = rand(N^2); W = W+W'; % random problem instance

params = [];

% run DS+
[XprojDsPlus, XdsPlus, lbDsPlus, ubDsPlus] = dsstar(W, params, 'ds+');

% run DS++
[XprojDsPlusPlus, XdsPlusPlus, lbDsPlusPlus, ubDsPlusPlus] = dsstar(W, params, 'ds++');

% run DS*
[XprojDsStar, XdsStar, lbDsStar, ubDsStar] = dsstar(W, params, 'ds*');

disp(['DS+ (lower bound, upper bound): ' num2str(lbDsPlus) ', ' num2str(ubDsPlus)]);
disp(['DS++ (lower bound, upper bound): ' num2str(lbDsPlusPlus) ', ' num2str(ubDsPlusPlus)]);
disp(['DS* (lower bound, upper bound): ' num2str(lbDsStar) ', ' num2str(ubDsStar)]);

