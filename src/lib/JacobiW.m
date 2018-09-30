function v1 = JacobiW(A,b,v0,n,w)
% Inputs:
%   A: HWxHW (typically sparse) matrix of the proble. 
%   b: HWx1 column vector containing the data. 
%   v0: HWx1 colun vector containing the initial guess for the solution. 
%   n: (scalar) number of Jacobi iterations performed. 
%   w: (scalar) weight. 
% -- 
% Implements n steps of the weighted Jacobi iterative method on the 
% problem Ax = b. (For w = 1 we have standard Jacobi method).  
% --

% (i) Decompose A in upper triangular, lower triangular and diagonal
L = tril(A,-1); U = triu(A,1); D = diag(diag(A));

% (ii) Define the iteration matrix and pre compute the data term multip.
R   = (1-w)*speye(size(A)) - w*(D\(L+U));
wDb = w*(D\b);

for i = 1:n, 
    v1 = R*v0 + wDb; % weighted Jacobi iteration step
    v0 = v1;         % update v0
end
end