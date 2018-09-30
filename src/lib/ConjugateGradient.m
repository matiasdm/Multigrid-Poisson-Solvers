function u1 = ConjugateGradient(A,b,u0,numIter) 
% [uf,R] = ConjugateGradient(A,b,u0,parameters) 
% Solve the problem Ax = b; using the conjugate gradient descend method.
% Inputs,
%   - A: NxN sparse or full matrix (semidefinite positive matrix) 
%   - b: Nx1 vector
%   - u0: seed (initial guess)
%   - numIter: number of iterations, 
%   - parameters: struct that may contain,
% Outputs, 
%   - uf: solution of equation Ax = b; 

% init. 
iter    = 0; 
r0      = b-A*u0(:);
p0      = r0; 

while iter<numIter,
    
    alpha = r0'*r0 / (p0'*A*p0); % alfa coeficient     
    % Update u_n+1,
    u1 = u0 + alpha*p0;        
    
    % Update r, alpha and p
    r1   = r0 - alpha*A*p0;
    beta = r1'*r1 / (r0'*r0);
    p0   = r1 + beta*p0;
    r0   = r1;
    
    % Update the number of iteration and diff. 
    iter = iter+1;
    u0   = u1;
    
end
end
