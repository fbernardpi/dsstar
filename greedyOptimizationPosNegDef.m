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
function [d1,d2] = greedyOptimizationPosNegDef(W,kerA,N,maxIter,beta, proxOperator)
    if ( ~exist('beta', 'var') || isempty(beta) )
        beta = 0.2;
    end
    
    tau = 4;
    
    eta = 0.1;
    proxOperator = @(x) x/(1 + tau*eta);
    
    if ( ~exist('ourEigs', 'file') )
        useOurEigs = 0;
        warning('File ourEigs.m not available. Patch the file eigs.m according to eigs_to_ourEigs.patch for better performance');     
    else
        useOurEigs = 1;
    end
    
    % gradient
    d1 = zeros(N,1);
    d2 = zeros(N,1);
    
    Wsubspace = kerA'*W*kerA;
    
    eigsTol = 1e-4;
    maxit = 100;
    for i=1:maxIter
        % This is still the asumptotically most expensive part
        currentMatPsd = kerA'* ...
            (kron(diag(d1), speye(N)) + kron(speye(N), diag(d2)))* ...
            kerA;
        
        T0 = Wsubspace - currentMatPsd;
        opts = struct('issym', 1, 'isreal', 1, 'tol', eigsTol, ...
            'disp', 0, 'maxit', maxit);
        
        T0 = (T0+T0')/2;
        [U_min,lambda_min] = eigs(T0, 1,'SA',opts);
        
        clear T0;
        
        
        Vplus = reshape((kerA*U_min).^2, [N,N]);
        
        T1 = Wsubspace + currentMatPsd;
        clear currentMatPsd;
        
        if ( beta ~= 0 )
            opts = struct('issym', 1, 'isreal', 1, 'tol', eigsTol, ...
                'disp', 0, 'maxit', maxit);%, 'v0', U_max);
            % here we use a modified version of eigs that terminates earlier in order
            % to avoid very slow processing times for non-simple spectra
            if ( useOurEigs )
                [U_max,lambda_max,~] = ourEigs((T1+T1')/2, 1,'LA',opts);
            else 
                [U_max,lambda_max,~] = eigs((T1+T1')/2, 1,'LA',opts);
            end
                
            if ( isnan(lambda_max) )
                lambda_max = 0;
                % warning('eigs() has not converged. Setting sigma_max = 0');
            end
            clear T1;
            
            Vminus = reshape((kerA*U_max).^2, [N,N]);
        else
            Vminus = 0;
            lambda_max = 0;
        end
        
        
        d1 = (d1+(1-beta)*tau*lambda_min*sum(Vplus,1)' ...
            -beta*tau*lambda_max*sum(Vminus,1)');
        d1 = proxOperator(d1);
        
        d2 = (d2+(1-beta)*tau*lambda_min*sum(Vplus,2) ...
            -beta*tau*lambda_max*sum(Vminus,2));
        d2 = proxOperator(d2);
        
        d1(d1>0) = 0;
        d2(d2>0) = 0;
    end