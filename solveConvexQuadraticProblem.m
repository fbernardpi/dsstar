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
function [Xsol,W,lambdaMinSubspace,z, H] = solveConvexQuadraticProblem(W, kerA, N, ...
        qpsolver,linTerm)
    
    if ( ~exist('qpsolver', 'var') )
        qpsolver = 'quadprog';
    end
    
    vec = @(x) x(:);
    
    % use QP
    lambda = vec(speye(N));
    
    x0 = lambda;
    
    H = kerA'*W*kerA;
    H = (H+H')/2;
    
    lambdaMinSubspace = eigs(H, 1,'sa')-1e-6;
    
    
    H = H - lambdaMinSubspace*speye(size(H));
    W = W - lambdaMinSubspace*speye(size(W));
    
    switch qpsolver
        case {'FW', 'fw'}
            Xsol = frankWolfeQuadratic(W, eye(N));
            z = [];
        case 'quadprog'
            if ~exist('linTerm', 'var')
                linTerm=zeros(size(W,1),1);
            end
            f = kerA'*W*lambda + 0.5*kerA'*linTerm;
            Aineq = -kerA;
            bineq = lambda;
            
            opts.Display = 'off';
            z = quadprog(H, f, Aineq, bineq, [], [], [], [], x0, opts);
            if ( isempty(z) )
                error('quadprog failed');
            end
            
            x = lambda + kerA*z;
            Xsol = reshape(x, N, N);
    end
end