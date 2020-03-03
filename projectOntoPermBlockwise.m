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
function Pproj = projectOntoPermBlockwise(Pin, D, nRows, nCols)
	if (~exist('nRows', 'var') )
		nRows = 1;
	end
	if (~exist('nCols', 'var') )
		nCols = 1;
	end
	Pin = sparse(Pin);
	Pproj = sparse(D*nRows, D*nCols);
	for c=1:nCols
		colIdx = (c-1)*D+1:c*D;
		for r=1:nRows
			rowIdx = (r-1)*D+1:r*D;
			
			currBlock = Pin(rowIdx,colIdx);
			currBlock = currBlock - min(currBlock(:)) + sparse(1);
			currBlock = currBlock*1e6;
			[~,currP] = sparseAssignmentProblemAuctionAlgorithm(currBlock);
			
			Pproj(rowIdx,colIdx) = currP;
		end
	end
end