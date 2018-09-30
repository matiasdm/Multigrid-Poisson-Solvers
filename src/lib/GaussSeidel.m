function v1 = GaussSeidel(A,b,v0,n)
% Inputs:
%   A: HWxHW (typically sparse) matrix of the problem. 
%   b: HWx1 column vector containing the data. 
%   v0: HWx1 colun vector containing the initial guess for the solution. 
%   n: (scalar) number of Jacobi iterations performed. 
% -- 
% Implements n steps of the weighted Jacobi iterative method on the 
% problem Ax = b. (For w = 1 we have standard Jacobi method).  
% --

% (i) Decompose A in upper triangular, lower triangular and diagonal
L = tril(A,-1); U = triu(A,1); D = diag(diag(A));

% (ii) GaussSeidel Iteration
%R   = -(D+L)\U; 
%DLb = (D+L)\b;
DL = (D+L);

for i = 1:n, 
    v1 = DL\(-U*v0 + b); 
    v0 = v1;         % update v0
end
end