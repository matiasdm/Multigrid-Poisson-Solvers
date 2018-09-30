function [v,flag,t] = PoissonInpainting(f,g,phi,v0,par)
% [v,flag,t] = PoissonInpainting(f,g,phi,v0,parameters) solves the problem:
%                 ----
%   lap(v(x,y)) = f(x,y)     on Omega 
%   v(x,y)      = g(x,y)     if (x,y) belongs to phi 
%   neumann BC               on dOmega.
%                 ----
% The problem is solved writting the discrete problem using a Self-adjoint 
% discrete formulation as explained in ref [1]. By default a Multi grid
% approach is used to solve the sparse linear problem, however GaiussSeidel
% or Conjugate gradient can also be selected. 
%
% Inputs: 
%   f: HxWxC Color or gray image with the laplacian data term 
%   g: HxWxC Color or gray image with dirichlet BC
%   phi: HxWxC binary image that corresponds to the indicator func of Phi. 
%   v0: [def matrix full of 0] initial guess of the solution.
%   parameters is a struct that may contain the following opt. params.:
%       .solver: [def 'MG'] can be: 'MG': multigrid, 'GS': Gauss-Seidel,
%                            'CG': conjugate gradients; 'BS': backslask
%       .verbose: [def 0] 1 display some text, 2 also some graphics  
%       .normRes: [def H*W] norm of the residual for which we considered we
%                 converged to the solution. 
%       .timeMax: [def 30] time in seconds after which iteration is
%                 interrupted. When this happend flag == 0, and we consider
%                 convergence was not reached.
%       .iter:    [def inf] break after these many iterations (not used by BS)
%       .numIter: [def 3] For GS and CG sets the number of iterations
%                 steps, for MG approaches, sets the number of pre- post-
%                 smoothing iterations. 
%       --- Additional params for 'MG',
%         .nu: [def 1] number of times a coarse layer is called (nu=1
%                      corresponds to a V-cycle, and nu=2 to a W-cycle)
%
% Output:
%   v: solution of the poisson problem defined above. 
%   flag: indicates if convergence has reached (flag = 1). 
%   t: time until convergence.
% ----------------------------------------------------------------------- %
% Ref:
% [1] MultiGrid Poisson Solvers, IPOL
% ----------------------------------------------------------------------- %
% matiasdm@fing.edu.uy                                        Paris 2017  %
% ----------------------------------------------------------------------- %

% Handle default parameter
if ~exist('par','var'),  % no optional parameters are provided 
   par = {};
end
if ~exist('v0','var'),   % no optional initialization is provided 
   v0 = g * 0.0;
end

% First check if the image is color or gray, if it is a color image, solve
% the problem for each color channel in a independen way:
[H,W,C] = size(f);
if C>1, % then we have color data,
    v  = zeros(H,W,C); % init
    flag = [0 0 0];
    t    = [0 0 0];
    for k = 1:C, % use parfor here for a par. computation. 
        disp(['Channel: ', num2str(k)])
        [v(:,:,k),flag(k),t(k)] = PoissonInpainting(...
                   f(:,:,k),g(:,:,k),phi(:,:,k),v0(:,:,k),par);
    end
    flag = min(flag); t = max(t);
    return % and quit the function
end
    
% From now on we are solving the problem assuming we have a gray image 

% -------------------------------------------- %
% -- Read optional input parameters         -- %
% -------------------------------------------- %
if isempty(par), % no optional parameters are provided, 
    nu      = 1; % set default values, 
    solver  = 'MG';
    verbose = 0;
    numIter = 3;
    timeMax = 30;
    normRes = H*W;
else
    % ==== solver ==== %
    if isfield(par,'solver'), solver = par.solver; else solver = 'MG'; end
    % ==== normRes ==== %
    if isfield(par,'normRes'), normRes = par.normRes; 
    else normRes = H*W; end
    % ==== timeMax ==== %
    if isfield(par,'timeMax'), timeMax = par.timeMax; 
    else timeMax = 30; end
    % ==== iter ==== %
    if isfield(par,'iter'), iter = par.iter; else iter = Inf; end
    % ==== nu ==== %
    if isfield(par,'nu'), nu = par.nu; else nu = 1; end
    % ==== numIter ==== %
    if isfield(par,'numIter'), numIter = par.numIter; else numIter = 3; end
    % ==== verbose ==== % 
    if isfield(par,'verbose'), verbose = par.verbose; else verbose = 0; end
end

%% (1) Build the lin sys. associated to the problem, 
% Compute the laplacian matrix.
HW = H*W; 
L  = createLaplacianOperator(H,W); 
P  = sparse(1:HW,1:HW,phi(:),HW,HW); % Sparse diag. matrix with bound. pos.
Pc = sparse(1:HW,1:HW,1-phi(:),HW,HW); % Sparse diag. matrix with bound. pos.
A  = -Pc*L*Pc + P; % Self adjoint-extended formulation (Eq. (15))
b  = Pc*(-f(:) + L*g(:));
v  = v0(:)-g(:); % init the solution as the seed, and apply the change of
                 % variable that renders the formulation self-adjoint

%% (2) Solve the lin sys Ax = b,
              
res  = normRes + 1; % init. 
tic; t = toc; % init.

switch solver, 
    case 'MG', % use Multi-Grid method,
        disp('running MG')
        parMG.verbose = verbose; 
        parMG.nu      = nu; 
        parMG.numIter = numIter;
        i=0;
        while res>normRes,
            v   = MG(A,b,v,H,W,parMG);
            res = norm(A*v-b);
            t   = toc;
            if t>=timeMax, disp('Timeout'); break; end;
            if i>=iter, break; end;
            i=i+1;
        end
        disp(['iter = ',num2str(i), '  t = ',num2str(t), '  res = ', num2str(res)])
    case 'GS', % use Gauss-Seidel method, 
        disp('running GS')
        i=0;
        while res>normRes,
            v   = GaussSeidel(A,b,v,numIter);
            res = norm(A*v-b);
            t   = toc;
            if t>=timeMax, disp('Timeout'); break; end;
            if i>=iter, break; end;
            i=i+1;
        end
        disp(['iter = ',num2str(i), '  t = ',num2str(t), '  res = ', num2str(res)])
    case 'CG', % use Conjugate gradient method,  
        disp('running CG')
        i=0;
        while res>normRes,
            v   = ConjugateGradient(A,b,v,numIter);
            res = norm(A*v-b);
            t   = toc;
            if t>=timeMax, disp('Timeout'); break; end;
            if i>=iter, break; end;
            i=i+1;
        end
        disp(['iter = ',num2str(i), '  t = ',num2str(t), '  res = ', num2str(res)])
    case 'BS', % use matlab backslash
        disp('running BS')
        v   = A\b; 
        res = norm(A*v-b);
        t   = toc;
        disp(['t = ',num2str(t), '  res = ', num2str(res)])
end
v    = reshape(v,[H W]); % vector --> matrix 
v = v + g;          % undo the change of variables used in Eq. (14-15)
flag = res<normRes; % flag == 1 if convergence was reached (red<normRes)

end
