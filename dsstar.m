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

function [Xsol, Xcont, lb, ub, Wtilde, lbOffset, lambdaMinSubspace, D1From, D1To, ...
        D2From, D2To, dFrom, dTo] = dsstar(W, params, initMode)
	N = sqrt(size(W,1));
    
    if ( ~exist('initMode', 'var') )
        initMode = 'ds*';
    end
    
    if ( ~exist('params', 'var') || ~isfield(params, 'qpsolver') )
        qpsolver = 'quadprog';
    else
        qpsolver = params.qpsolver;
    end
    
    if ( ~exist('params', 'var') || ~isfield(params, 'maxIterGreedyOpt') )
        maxIterGreedyOpt = 10;
    else
        maxIterGreedyOpt = params.maxIterGreedyOpt;
    end
    
    if ( exist('params', 'var') && isfield(params, 'useOrthogonalKerA') )
        useOrthogonalKerA = params.useOrthogonalKerA;
    else
        useOrthogonalKerA = 1;
    end
    
    if ( exist('params', 'var') && isfield(params, 'nStepsPf') )
        nSteps = params.nStepsPf;
    else
        nSteps = 10;
    end
    if ( exist('params', 'var') && isfield(params, 'stepMode') )
        stepMode = params.stepMode;
    else
        stepMode = 'lin';
    end
    if ( exist('params', 'var') && isfield(params,'useLambdaMinMaxInPf') )
        useLambdaMinMaxInPf = params.useLambdaMinMaxInPf;
    else
        useLambdaMinMaxInPf = 0;
    end

    switch stepMode
        case 'rnd'
            stepDeltas = cumsum(rand(nSteps,1));
            stepDeltas = stepDeltas./max(stepDeltas);
        case 'lin'
            stepDeltas = linspace(0,1,nSteps+1)';
            stepDeltas(1) = [];
    end
    
            
    fwFcnHandle = @frankWolfeQuadratic;
    
    
    
    
    A = [kron(speye(N), ones(1,N)); ...
    kron(ones(1,N), speye(N))];
    A(end,:) = [];

    switch initMode
        case 'ds+'
            b = ones(2*N-1,1);
            
            [~,lambdaMin] = eigs(W,1, 'sa');
            [~,lambdaMax] = eigs(W,1, 'la');
            Wtilde = W - lambdaMin*speye(N^2) + 1e-8;
            
            options.Display = 'off';
            z = quadprog(Wtilde, zeros(N^2,1), [], [], ...
                A, b, zeros(N^2,1), ones(N^2,1), [], options);
            Xcont = reshape(z, N, N);
            
            D1From = zeros(N,N);
            D2From = zeros(N,N);
            dFrom = lambdaMin*ones(N^2,1);
            
            D1To = zeros(N,N);
            D2To = zeros(N,N);
            dTo = lambdaMax*ones(N^2,1);
            
            lbOffset = N*lambdaMin;
            lambdaMinSubspace = nan;
        case 'ds++'
            kerA_ortho = null(full(A));

            [Xcont,Wtilde,lambdaMinSubspace,~,H] = ...
                solveConvexQuadraticProblem(W,kerA_ortho,N, qpsolver);
            lambdaMaxSubspace = eigs(H, 1, 'la');
            
            D1From = zeros(N,N);
            D2From = zeros(N,N);
            dFrom = lambdaMinSubspace*ones(N^2,1);
            
            D1To = zeros(N,N);
            D2To = zeros(N,N);
            dTo = lambdaMaxSubspace*ones(N^2,1);
            
            lbOffset = N*lambdaMinSubspace;
        case 'ds*'
            kerA_ortho = null(full(A));
            if ( useOrthogonalKerA )
                kerA = kerA_ortho;
            else % slightly faster
                kerA = nullDs(N, 0);
            end
            
            alpha = [];
            [D1From,D2From] = greedyOptimizationPosNegDef(...
                W, kerA, N, maxIterGreedyOpt, alpha);
            
            D1From = diag(D1From);
            D2From = diag(D2From);
            D1To = -D1From;
            D2To = -D2From;
            
            % construct final matrix
            Wtilde = W - kron((D1From), speye(N)) - kron(speye(N), (D2From));
            
            [Xcont,Wtilde,lambdaMinSubspace] = ...
                solveConvexQuadraticProblem(Wtilde,kerA_ortho,N,qpsolver);
            
            dFrom = zeros(N^2,1);
            dTo = dFrom;

            if ( useLambdaMinMaxInPf )
                Hpsd = kerA'*(W - kron((D1From), speye(N)) - kron(speye(N), (D2From)))*kerA;
                Hnsd = kerA'*(W - kron((D1To), speye(N)) - kron(speye(N), (D2To)))*kerA;
                Hpsd = (Hpsd+Hpsd')/2;
                Hnsd = (Hnsd+Hnsd')/2;
                lambdaMinSubspace = eigs(Hpsd, 1, 'sa');
                lambdaMaxSubspace = eigs(Hnsd, 1, 'la');
                dFrom = lambdaMinSubspace*ones(N^2,1);
                dTo = lambdaMaxSubspace*ones(N^2,1);
            end
            
            clear kerA_ortho;
            lbOffset = trace(D1From)+trace(D2From) + lambdaMinSubspace*N;
    end


    currSol = Xcont;

    % run path following
    for i=1:nSteps
        currAlpha = stepDeltas(i);
        
        % compute current W
        currD1 = (1-currAlpha)*D1From + currAlpha*(D1To);
        currD2 = (1-currAlpha)*D2From + currAlpha*(D2To);
        currD = (1-currAlpha)*dFrom + currAlpha*(dTo);
        
        currW = W - kron(currD1, speye(N)) - ...
            kron(speye(N), currD2) - diag(currD);
        currW = (currW+currW')/2;

        % local optimisation 
        currSol = fwFcnHandle(currW, currSol, currD);
        
        % convergence check
        currSolProj = projectOntoPermBlockwise(currSol, N);
        if ( abs(trace(currSol*currSolProj') - N) < 1e-4 )
            break;
        end
    end
    
    Xsol = currSolProj;
    
    lb = Xcont(:)'*Wtilde*Xcont(:) + lbOffset;
    ub = Xsol(:)'*W*Xsol(:);
end