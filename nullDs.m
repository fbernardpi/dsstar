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

function kerDs = nullDs(N)
    % here we construct a sparse matrix kerDs that spans the nullspace of A,
    % where A is the linear equality constraint matrix such that A*vec(X)
    % encodes that the rows and columns of X sum to 1
    A = [kron(speye(N), ones(1,N)); ...
        kron(ones(1,N), speye(N))];
    A(end,:) = []; % remove linear dependency
    
    pairs = [(1:N-1)' (2:N)'];
    
    dimKer = N^2 - (2*N-1);
    kerDs = sparse(N^2, dimKer);
    colIdx = 1;
    for i=1:size(pairs,1)
        for j=1:size(pairs,1)
            % TODO: implement non-forloop vectorised version
            pairI = pairs(i,:);
            pairJ = pairs(j,:);
            
            xi = sparse(pairI, [1 1], [1 -1], N, 1);
            yi = sparse(pairJ, [1 1], [1 -1], N, 1);
            
            kerDs(:,colIdx) = kron(xi,yi)./2; % divide by 2 to obtain unit norm
            colIdx = colIdx + 1;
        end
    end

    assert(norm(A*kerDs, 'fro') < 1e-6);
end