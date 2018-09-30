% ============================================================ %
% MultiGrid poisson solvers, 
% ------------------------------------------------------------ %
% This code illustrates who to use the codes
% and techiques described in the artilcle:
%
% [Di Martino et al. 2017]
% M. Di Martino and G. Facciolo 
% "Multigrid Poisson Solvers", Image Processing On Line IPOL, 2017.
% ------------------------------------------------------------ %
% copyright (c) 2017,
% Matias Di Martino <matiasdm@fing.edu.uy>
% Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
%
% Licence: This code is released under the AGPL version 3.
% Please see file LICENSE.txt for details.
% ------------------------------------------------------------ %
% Comments and sugestions are welcome at: matiasdm@fing.edu.uy
% M. Di Martino and G. Facciolo 
% ============================================================ %

close all
clear all
clc
home

% add the path where the important functions are and
% also the folder that contains some example images.
addpath src src/lib

% &&  Menu  &&  
% Example 1: we will compare solvers. 
% Example 2: we will compare different formulations of the discrete problem.   
% &&&&&&&&&&&& 

% ============================================== %
%%  Example 1: Interpolation from edges         %%
% ============================================== %
% In this first illustrative example, we will focus on comparing different
% iterative solvers of the poisson equations. 
fprintf(' ======================================= \n')
fprintf(' = Example 1: interpolation from edges = \n')
fprintf(' ======================================= \n')

% Read a test input image, 
filename = 'images/baboon.png';
I        = double(imread(filename)); 
[H,W,C]  = size(I); 

figure(1), imshow(uint8(I));    
% Shortcut to partial derivative, 
dx = @(U) [U(:,2:end,:) - U(:,1:end-1,:)  0*U(:,1,:)]; % Forward fin diff. 
dy = @(U) [U(2:end,:,:) - U(1:end-1,:,:); 0*U(1,:,:)]; % with Neumman BC.

% Keep some edges and interpolate the rest of the image with lap. Eq. 
th   = 90; % <-- #### MODIFY THIS VALUE TO CHANGE THE DOMAIN PHI #######
Mask = ( dx(I).^2 + dy(I).^2 ) > th^2;
Mask = max(Mask,[],3);
Mask = cat(3,Mask,Mask,Mask);

% ------------------------------------------------------ %
% Define the poisson inpainting problem as:              %
%   lap(v(x,y)) = f(x,y)     on Omega                    %
%   v(x,y)      = g(x,y)     if (x,y) belongs to phi     %
%   neumann BC               on dOmega.                  %
phi  = Mask;                                             %
f    = 0*ones(H,W,C);                                    %
g    = zeros(H,W,C); g(Mask) = I(Mask);                  %
u    = I;                                                %
v0   = zeros(H,W,C);                                     %
% ------------------------------------------------------ %

% --- display problem input --- %
figure(2),
subplot(131), imshow(uint8(255*phi)), title('\Phi'), axis image, 
subplot(132), imshow(uint8(f)),   title('f'), axis image, 
subplot(133), imshow(uint8(g)),   title('g'), axis image

% ============================================== %
% Test diffrent iterative Solvers                %
% ============================================== %
poissonPar.normRes = H*W*1e-7; % min norm for convergence
poissonPar.timeMax = 5; % max. time one have to reach convergence
poissonPar.numIter = 3;

% ----------------------------------------------------------- % 
% (C.1) GaussSeidel Method 
% ----------------------------------------------------------- %
poissonPar.solver = 'GS';
[v_gs,f_gs,t_gs]  = PoissonInpainting(f,g,phi,v0,poissonPar);
e_gs = norm(v_gs(:)-u(:)); figure, imshow(uint8(v_gs)), title('GS')

% ----------------------------------------------------------- % 
% (C.2) Conj gradient
% ----------------------------------------------------------- %
poissonPar.solver = 'CG';
[v_cg,f_cg,t_cg]  = PoissonInpainting(f,g,phi,v0,poissonPar);
e_cg = norm(v_cg(:)-u(:)); figure, imshow(uint8(v_cg)), title('CG')

% ----------------------------------------------------------- % 
% (C.3) Multi grid method
% ----------------------------------------------------------- %
poissonPar.solver = 'MG';
poissonPar.nu     = 1; % 1 = V-scheme, 2 = W-scheme
[v_mg,f_mg,t_mg] = PoissonInpainting(f,g,phi,v0,poissonPar);
e_mg = norm(v_mg(:)-u(:)); figure, imshow(uint8(v_mg)), title('MG')

% ============================================== %
%     Display results                            %
% ============================================== %
% See which was the fastest method 
mintime = 1;%min([t_gs(end); t_cg(end); t_mg(end)]);
fprintf(' =================================================== \n')
fprintf('    + GaussSeidel : error = %5.1e , time %3.1f, converged = %d\n', ...
         norm(v_gs(:)-u(:))/numel(u), t_gs(end)/mintime, f_gs )
fprintf('    + Conj.Gradi. : error = %5.1e , time %3.1f, converged = %d\n', ...
         norm(v_cg(:)-u(:))/numel(u), t_cg(end)/mintime, f_cg )
fprintf('    + Multi Grid  : error = %5.1e , time %3.1f, converged = %d\n', ...
         norm(v_mg(:)-u(:))/numel(u), t_mg(end)/mintime, f_mg )
fprintf(' =================================================== \n')
%  ------------------------------------- % 
clear all 

% ================================================ %
%%  Example 2: Image blending in the grad, domain %%
% ================================================ %
fprintf(' =================================================== \n')
fprintf(' = Example 2: Contrast enhancement/reduction       = \n')
fprintf(' =================================================== \n')
I_original = double(imread('images/catedral.png'));
hsv_im = rgb2hsv(I_original);
I = hsv_im(:,:,3); % Enhance the value and keep original colors.
[H,W,~] = size(I); 
Phi = I>40; % We will define the domain Phi as the set of bright pixels

% Compute the laplacian of the input image, 
L = createLaplacianOperator(H,W);
lap_I = reshape(L*I(:),[H W]); 

% To illustrate a practical application that may be solved using the
% poisson equation, we will illustrat how the dark portions of the image
% can be enhanced by amplifing the laplacian of the image in the interior
% of the dark areas, while keeping the original values of the image in the
% bright areas. To that end, f (the laplacian term) is defined as the
% original laplacian lap_I multiplied by a scalar factor alpha. g is
% set to the original values of the image in phi. 
alpha = 1; % amplification factor. 
f = zeros(H,W); g = zeros(H,W); % init. 
% ------------------------------------------------------ %
% Define the poisson inpainting problem as:              %
%   lap(v(x,y)) = f(x,y)     on Omega                    %
%   v(x,y)      = g(x,y)     if (x,y) belongs to phi     %
%   neumann BC               on dOmega.                  %
g(Phi)  = I(Phi);                                        %
f(~Phi) = alpha*lap_I(~Phi);                             %
v0      = zeros(H,W);                                  %
% ------------------------------------------------------ %

% --- display problem input --- %
figure(1),
subplot(131), imshow(uint8(255*Phi)), title('\Phi'), axis image, 
subplot(132), imshow(f,[]),   title('f'), axis image, 
subplot(133), imshow(uint8(g)),   title('g'), axis image

% -------------------------------------------------------- %
% (1) Build the Self-Adjoint Extended Formulation 
% -------------------------------------------------------- %
HW = H*W;
P  = sparse(1:HW,1:HW,Phi(:),HW,HW); % Sparse diag. matrix with bound. pos.
Pc = sparse(1:HW,1:HW,1-Phi(:),HW,HW);
A  = -Pc*L*Pc + P; % Self adjoint-extended formulation (Eq. (16))
b  = Pc*(-f(:) + L*g(:)) + P*g(:);  % **
% Note: the data term here already include the correction of the change of
% variable u' = u-g. 

v  = v0(:); % init the solution as the seed, 

parMG.verbose = 0; 
parMG.nu      = 1; 
parMG.numIter = 1;
% Solve the sparse lin sys. 
%v = MG(A,b,v,H,W,parMG);  % <-- using the multi-grid method presented 
v = A\b; % <-- using matlab backslash


% Normalization 
v = 255 * ( v-min(v(:)) ) / ( max(v(:)) - min(v(:)) );
Isa_ext = hsv2rgb(cat(3,hsv_im(:,:,1),hsv_im(:,:,2),reshape(v,[H W])));
imwrite(uint8(Isa_ext),'Isa_ext.png')

% -------------------------------------------------------- %
% (2) Build the Extended Formulation 
% -------------------------------------------------------- %
% -- Multi - Grid solver -- %
A  = -Pc*L + P; % Self adjoint-extended formulation (Eq. (16))
b  = Pc*(-f(:)) + P*g(:);

v  = v0(:); % init the solution as the seed, 

% Solve the sparse lin sys. 
%v = MG(A,b,v,H,W,parMG);  % <-- using the multi-grid method presented 
v = A\b; % <-- using matlab backslash


% Normalization, 
v = 255 * ( v-min(v(:)) ) / ( max(v(:)) - min(v(:)) );
Iext = hsv2rgb(cat(3,hsv_im(:,:,1),hsv_im(:,:,2),reshape(v,[H W])));
imwrite(uint8(Iext),'Iext.png')

figure(2), 
subplot(131), imshow(uint8(I_original)), title('input image'); 
subplot(132), imshow(uint8(Isa_ext)), title('Isa_ext'); 
subplot(133), imshow(uint8(Iext)), title('Iext'); 

% -- Sanity Check -- % 
% You can set alpha = 1, and verify if the original image is recover. Try
% this using backslash and multigrid method, and see if both formulation of
% the problem (mathematically equivalent) are equally suitable for
% multigrid iterative solvers... 
fprintf('<|I_original - I_sa-ext|> = %3.1f \n ', ...
    mean(abs(double(I_original(:))-double(Isa_ext(:)))));
fprintf('<|I_original - I_ext|> = %3.1f \n ', ...
    mean(abs(double(I_original(:))-double(Iext(:)))));
% ------------------ %
