function K = kron(A,B)
% Overload default MATLAB implementation to support OptiSpline Function
% objects.

% This file is part of gssshinfcd.
% Copyright (c) 2024, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% gssshinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% gssshinfcd is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
% License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with gssshinfcd. If not, see <https://www.gnu.org/licenses/>.

       cA = mat2cell(A,ones(size(A,1),1),ones(size(A,2),1));
       K = cellfun(@(x) x*B, cA, 'un', 0);
       K = reshape(vertcat(K{:}),[size(A,1)*size(B,1),size(A,2)*size(B,2)]);
   
end