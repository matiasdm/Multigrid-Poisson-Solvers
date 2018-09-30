function runMG(finput, fmask, foutput, nu, smooth_iter, iter, solver, lap_multiplier, timeout)

if ~exist('finput','var'), finput='input_0.png'; warning('Parameter finput missing (automatically set to input_0.png)'); end;
if ~exist('fmask','var'), fmask='input_1.png'; warning('Parameter finput missing (automatically set to input_1.png)'); end;
if ~exist('foutput','var'), foutput='output_0.png'; warning('Parameter foutput missing (automatically set to output_0.png)'); end;
if ~exist('nu','var'), nu=1; warning('Parameter nu missing (automatically set to 1)'); else nu=str2num(nu); end;
if ~exist('smooth_iter','var'), smooth_iter=3; warning('Parameter smooth_iter missing (automatically set to 3)'); else smooth_iter=str2num(smooth_iter); end;
if ~exist('iter','var'), iter=5; warning('Parameter iter missing (automatically set to 5)'); else iter=str2num(iter); end;
if ~exist('solver','var'), solver='MG'; warning('Parameter solver missing (automatically set to MG)'); end;
if ~exist('lap_multiplier','var'), lap_multiplier=2; warning('Parameter lap_multiplier missing (automatically set to 2)'); else lap_multiplier=str2num(lap_multiplier); end;
if ~exist('timeout','var'), timeout=10; warning('Parameter timeout missing (automatically set to 10s)'); else timeout=str2num(timeout); end;

if exist(finput, 'file')
    % Matlab is able to find out type without extension
    try g = imread(finput);
    catch
        % input_0 is not an image
        error(['Unable to read file: ' finput]);
    end
    
end
if exist(fmask, 'file')
    % Matlab is able to find out type without extension
    try phi = imread(fmask);
    catch
        % input_1 is not an image
        error(['Unable to read file: ' fmask]);
    end
    if(size(phi,1)~=size(g,1) | size(phi,2)~=size(g,2) )
        error(['Mask is not the same size as the image.']);
    end
    
end

% normalize channels
if(size(g,3) ~= size(phi,3))
   phi=repmat(phi(:,:,1),[1,1,size(g,3)]);
end

g=double(g);
phi=double(255-phi); % invert the mask


% compute the inputs
f=g*0;
L= createLaplacianOperator(size(g,1),size(g,2));
for i =1:size(g,3),
   tmp = g(:,:,i);
   tmp = L*tmp(:);
   f(:,:,i) = reshape(tmp,[size(g,1), size(g,2)]);
end
f = f*lap_multiplier;
v0=g*0;

% parameters 
par.iter = iter; % fix the iterations
par.normRes = 0; % do not stop because of this
par.timeMax = timeout;% usually 10 seconds per channel < 30 (ipol usability limit)
par.solver = solver; 
par.numIter = smooth_iter;
par.nu = nu;
% run
[v,flag,t] = PoissonImpainting(f,g,phi,v0,par);



imwrite(uint8(v), [foutput]);
imwrite(uint8(g), ['output_1.png']);
imwrite(uint8(phi), ['output_2.png']);

% remove MG cache from the server
if exist('precompRP','dir'), rmdir('precompRP','s'); end

end



