function L = createLaplacianOperator(H,W) 
% L = createLaplacianOperator(H,W)
% This function creates the 5 points stencil laplatian operator ASSUMING
% NEUMANN BOUNDARY CONDITIONS. H and W corresponds to the size of the first
% and second dimension of the problem. 

Mask = ones(H,W);
Mask = padarray(Mask,[1 1], 0, 'both');
idx  = find(Mask(:)); % Index of interior points

[Hp,Wp] = size(Mask); % size of the padded domain.
N       = Hp*Wp;      % num of px. (in the ext. domain). 

L =   sparse(idx,idx   ,-4,N,N) ... % :central px {-4*I[i,j] = -4I(x,y)}
    + sparse(idx,idx+1 , 1,N,N) ... % :frw y-par der. {I[i+1,j] = I(x,y+1)}
    + sparse(idx,idx-1 , 1,N,N) ... % :bkw y-par der. {I[i-1,j] = I(x,y-1)}
    + sparse(idx,idx+Hp, 1,N,N) ... % :frw x-par der. {I[i,j+1] = I(x+1,y)}
    + sparse(idx,idx-Hp, 1,N,N);    % :bkw x-par der. {I[i,j-1] = I(x-1,y)}

% Keep just the equation that correspond to points in the domain, 
L = L(idx,idx);

clear Hp Wp N Mask idx 
N = H*W; % number of pixels (without padding). 

% Correct the weight of px that belong to the border, 
L = L - sparse(1:N, 1:N, sum(L,2), N, N);
end