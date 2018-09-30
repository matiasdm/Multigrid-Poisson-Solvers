% /////////////////////////////////////////////////////////////////////// %
% // Mulrigrid solver                                                  // %
% /////////////////////////////////////////////////////////////////////// %
function v1 = MG(A,b,v0,H,W,par)
% Inputs: 
%   A: HWxHW (typically sparse) matrix of the problem. 
%   b: HWx1 column vector containing the data term. 
%   v0: HWx1 colun vector containing the initial guess for the solution. 
%   par: a struct that may contain additional parameters, for example:
%       .nu: [def 1] number of times a coarse layer is called (nu=1
%            corresponds to a V-cycle, and nu=2 to a W-cycle)
%       .numIter:    [def 3] num of pre smoothing relaxation 
%       .numIterPos: [numIter] num of post smoothing, by defaul is equal 
%                    to the number of pre smoothing. 
%       .verbose: [def 0] 1- displays some text, 2- also some graphics,%
%
% -- 
% This function implements MG method described in the manuscript.   
% --

% -------------------------------------------- %
% -- Read optional input parameters         -- %
% -------------------------------------------- %
if isempty(par), % no optional parameters are provided, 
    nu         = 1; % set default values, 
    verbose    = 0;
    numIter    = 3;
    numIterPos = numIter;
else
    % ==== nu ==== %
    if isfield(par,'nu'), nu = par.nu; else nu = 1; end
    % ==== numIter ==== %
    if isfield(par,'numIter'), numIter = par.numIter; else numIter = 3; end
    % ==== numIter ==== %
    if isfield(par,'numIterPos'), numIterPos = par.numIterPos; 
        else numIterPos = numIter; end
    % ==== verbose ==== % 
    if isfield(par,'verbose'), verbose = par.verbose;
    else verbose = 0; end
end


% As MG multi-grid solver is called in a iterative (nested) way, we need to
% treat two different cases separatelly, (1) when we reached the bottom
% layer (e.g. the size of the discrete domain is H = 2, W = 2) then we just
% solve the problem, (2) otherwise we continue the iterative "two-grid"
% scheme process. 
if max([H; W])<=2, % then solve the problem 
    v1 = A\b; % any num solver can be used here,
    if verbose>0, fprintf('E'); end
else % perform a two-grid correction scheme, 

% Padd inputs if its size is not multiple of 2,
[Id_padding,Id_unpadding] = createPaddingMatrices(H,W);
% if H and W are multiple of 2, then padding and unpadding images are the
% identity (and do not have any effect). 

% Make the domain dimension multiple of 2
if mod(H,2)~=0, H = H+1; end, if mod(W,2)~=0, W = W+1; end, 

% First compute restriction and interpolation ops. 
[I2h_h,Ih_2h] = DefineInterpolationAndRestrictionOp(H,W);

% Combine interpolation & restriction with padding & unpadding
I2h_h = I2h_h*Id_padding;
Ih_2h = Id_unpadding*Ih_2h;

% -------------------------------------------------- %
% -- (1) n1 relaxation steps in the fine grid     -- %
% -------------------------------------------------- %
n1  = numIter; % pre-smoothing iteration steps 
if n1>0, vh  = GaussSeidel(A,b,v0,n1); else vh = v0; end
if verbose>0, fprintf('.'), end

% -------------------------------------------------- %
% -- (2) compute the residual after relaxation    -- %
% -------------------------------------------------- %
rh  = b - A*vh; 

% -------------------------------------------------- %
% -- (3) Restrict the residual to a coarser layer -- %
% -------------------------------------------------- %
r2h = I2h_h*rh;

% -------------------------------------------------- %
% -- (4) Relax on the coarser grid                -- %
% -------------------------------------------------- %
% In the two-grid approach, here we just perform some relaxation steps,
% i.e,
% e2h = GaussSeidel(I2h_h*A*Ih_2h,r2h,nIter);
%                   = A2h       = b2h = numberOfRelaxationSteps

% However, it is more efficient to iteratively implement the two-grid correction
% scheme to solve the "new" problem Ax = b where, A = A2h, and b = r2h.
A2h = I2h_h*A*Ih_2h;
b2h = r2h;
e2h  = 0*r2h; % zero initial guess,
for k = 1:nu, 
    if verbose>0, fprintf('\\'), end
    e2h = MG(A2h, b2h, e2h, H/2, W/2 ,par); % nu call of the coarse layer
%   v1  = solver(A, b, v0, H, W, params);   %(nu = 1//2 => V//W-cycle)   
    if verbose>0, fprintf('/'), end
end

% -------------------------------------------------- %
% -- (5) Interp. the sol. back to the fine grid   -- %
% -------------------------------------------------- %
eh  = Ih_2h*e2h;

% -------------------------------------------------- %
% -- (6) n2 post-smoothing relaxation steps       -- %
% -------------------------------------------------- %
vh = vh + eh; % up. the sol. with the cor. factor obt. in the coarse layer  
n2 = numIterPos; % number of post-smoothing steps, 
if n2>0, v1 = GaussSeidel(A,b,vh,n2); else v1 = vh; end
if verbose>0, fprintf('.'), end

end %ELSE

end

function [Id_padding,Id_unpadding] = createPaddingMatrices(H,W),
%% padd and unpadd an additional column/row (as a matrix multiplication) 
%% so as to be divisible by two

Wo = W; Ho = H; % Input size,
sho = [Ho,Wo];
HWo = prod(sho);

if mod(W,2)~=0, % padd adding a column,
       W = W+1;
end
if mod(H,2)~=0, % padd adding a row,
           H = H+1;
end
sh = [H,W]; % output size 
HW = prod(sh);

% ====================================== %
% Padding matrix                         %
% ====================================== %
[py,px] = ind2sub( sh, [1:HW]);           % get pixel indices for the output image
i = sub2ind(sh , py, px);                 % get rows indices of the operator (y,x) (identity)
j = sub2ind(sho, min(py,Ho), min(px,Wo)); % get corresponding input columns (min(y,Ho), min(x,Wo)) 
Id_padding   = sparse(i,j,1,HW,HWo);

% ====================================== %
% UnPadding matrix                       %
% ====================================== %
[py,px] = ind2sub(sho , [1:HWo]);         % get pixel indices for the output image
i = sub2ind(sho, py, px);                 % get rows indices of the operator (y,x) (identity)
j = sub2ind(sh , py, px);                 % get corresponding input columns (y,x) but with the right shape
Id_unpadding = sparse(i,j,1,HWo,HW);
end


function [I2h_h,Ih_2h] = DefineInterpolationAndRestrictionOp(H,W)
% Define interpolation and restriction operators acording to the
% description provided in Sec 3.2.1 and algorithms 3 and 4. 
if exist(['precompRP/RP_' num2str(H) '_' num2str(W) '.mat'],'file'),
    % then load precomputed matrices to save time, 
    load(['precompRP/RP_' num2str(H) '_' num2str(W) '.mat']);
    % We read the restriction operator I2h_h and compute the interpolation
    % operator as 4 times the transpose of the restriction one. The
    % multiplicative factor is to mantain functions norm when operators are
    % applied. 
    Ih_2h = 4*I2h_h';            % Interpolation (2h --> h)
else % compute matrices,     
    HW = H*W;
% --- Definition of Basic Restriction Operator R0 --- % 
mask = zeros(H,W); % Image full of zeros, 
mask(1:2:end,1:2:end) = 1; % [1 0 1 0; 0 0 0 0; 1 0 1 0; 0 0 0 0]
ind = find(mask(:)); % pos where mask == 1;
M = numel(ind);
R1 = sparse(1:M,ind,1,M,HW);

mask = zeros(H,W); % Image full of zeros, 
mask(2:2:end,1:2:end) = 1; % [1 0 1 0; 0 0 0 0; 1 0 1 0; 0 0 0 0]
ind = find(mask(:)); % pos where mask == 1;
M = numel(ind);
R2 = sparse(1:M,ind,1,M,HW);

mask = zeros(H,W); % Image full of zeros, 
mask(1:2:end,2:2:end) = 1; % [1 0 1 0; 0 0 0 0; 1 0 1 0; 0 0 0 0]
ind = find(mask(:)); % pos where mask == 1;
M = numel(ind);
R3 = sparse(1:M,ind,1,M,HW);

mask = zeros(H,W); % Image full of zeros, 
mask(2:2:end,2:2:end) = 1; % [1 0 1 0; 0 0 0 0; 1 0 1 0; 0 0 0 0]
ind = find(mask(:)); % pos where mask == 1;
M = numel(ind);
R4 = sparse(1:M,ind,1,M,HW);

% ------------------------------------------ %
% --- Def of the smoothing kernel matrix --- %
% ------------------------------------------ %
% IMPOSING NEUMAN BOUNDARY CONDITIONS.       %
% ------------------------------------------ %

% prepare sparse matrix indices
sh = [H,W];
[py,px] = ind2sub( sh, [1:HW]);           % get pixel indices for the output image
i = sub2ind(sh , py, px);                 % get rows indices of the operator (y,x) (identity)

G = sparse(HW,HW);
% generate the offsets of the stencil
alldxdy = [[0,0,1];[1,0,1/2];[-1,0,1/2];[0, 1,1/2];[ 0,-1,1/2];
           [1,1,1/4];[-1,1,1/4];[1,-1,1/4];[-1,-1,1/4]];
for t =1:9,
   dy=alldxdy(t,1); dx=alldxdy(t,2); val=alldxdy(t,3);
   % get corresponding input columns mirroring at the boundaries ie max(min(ycoord,H),1) 
   j = sub2ind(sh, max(min(py+dy,H),1), max(min(px+dx,W),1)); 
   G = G + sparse(i,j,val,HW,HW);
end
G = G/4;

% ------------------------------------------ %
% --- Def of the smoothing kernel matrix --- %
% ------------------------------------------ %
% IMPOSING NEUMAN BOUNDARY CONDITIONS.       %
% ------------------------------------------ %
I2h_h = 1/4*(R1+R2+R3+R4)*G; % Restriction (h --> 2h)
Ih_2h = 4*I2h_h';            % Interpolation (2h --> h)

% ~~~ Sanity check, ~~~ % 
% One can verify that Ih_2h*I2h_h*C(:) = C(:) for any constant matrix 
% C = c*ones(H,W)
% ~~~~~~~~~~~~~~~~~~~~~ %

% Save them so we do not need to recompute them every time, 
if ~exist('precompRP','dir'); mkdir('precompRP'); end 
save(['precompRP/RP_' num2str(H) '_' num2str(W) '.mat'],'I2h_h');

end

end

% /////////////////////////////////////////////////////////////////////// %
% // Basic (one-grid) iterative solvers                                // %
% /////////////////////////////////////////////////////////////////////// %
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


% /////////////////////////////////////////////////////////////////////// %
% // Aux functions and definitions                                     // %
% /////////////////////////////////////////////////////////////////////// %
function normIh = norm_h(Ih,h)
% Implementation of the discrete L2 norm. For a mathematical definition of
% the norme implemented here, see [3,pg. 55].
% Ih denotes a two dimensional array to whom the norm is computed, and h
% accounts for the distance between pixels. 
%
normIh = sqrt( h^2 * sum( Ih(:).^2 ) );
end
