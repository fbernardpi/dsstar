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
function [X, objVals] = frankWolfeQuadratic(W, X0, d, convergenceThreshold, nMaxFWIters)

    N = size(X0,1);
    if ( ~exist('d', 'var') || isempty(d) )
        d = zeros(N^2,1);
    end
    
    if ( ~exist('convergenceThreshold', 'var') )
        convergenceThreshold = 1e-2;
    end
    if ( ~exist('nMaxFWIters', 'var') )
        nMaxFWIters = 1000;
    end
%     obj = @(X) X(:)'*W*X(:) + d(:)'*X(:);
    
    
    Xold = X0;

    objVals = nan;
    for jj=1:nMaxFWIters
        % find direction D
        % we want to solve D = argmax_{D =[D_1; ...; D_k] : D_i \in Perm}
        % < gradF(Pold), D>

        gradF = 2*W*Xold(:) + d;

        gradFScaled = -gradF;
        gradFScaled = gradFScaled - min(gradFScaled(:));
        gradFScaled = gradFScaled./max(gradFScaled(:));
        gradFScaled = gradFScaled*1e6;

        % we need to add 1 to make sure that there exists a feasible solution
        gradFScaled = reshape(gradFScaled, N, N)+1; 
        [~,D] = sparseAssignmentProblemAuctionAlgorithm(...
            sparse(gradFScaled), [], [], 0);

        DminusXold = D-Xold;

        % solve for step size (closed-form solution)
        a = trace(DminusXold(:)'*W*DminusXold(:));
        b = 2*trace(DminusXold(:)'*W*Xold(:)) + DminusXold(:)'*d;

        if ( abs(a) < 1e-8 )
            X = Xold;
            break;
        end
        bByMinus2a = -b/(2*a);

        if ( a > 0 )
            if ( bByMinus2a <= 0 )
                eta = 0;
            elseif ( bByMinus2a <= 1 )
                eta = bByMinus2a;
            else
                eta = 1;
            end
        else
            if ( bByMinus2a <= 0.5 )
                eta = 1;
            else
                eta = 0;
            end
        end


        X = Xold + eta*DminusXold;

        objVals(jj+1) = X(:)'*W*X(:);


        if ( abs(objVals(jj+1)/objVals(jj) - 1) < convergenceThreshold || ...
                abs(a) < 1e-8 )
            break;
        end
        Xold = X;
    end
end