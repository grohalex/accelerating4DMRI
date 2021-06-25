function [img,c,cost,imgsv,csv,p] = TrACR(data,p)

% correct for k-space trajectory errors and reconstruct an image
%
% data = k-space data (number of points x number of coils)
% p = parameter structure
%
% Copyright 2015, Julianna Ianni and Will Grissom, Vanderbilt University

% Do some argument processing
%if ~isfield(p,'useparfor')
%    p.useparfor = 0;
%end
if ~isfield(p,'bls')
    p.bls = 1;
end

% initialize k-space error guess
if isfield(p,'cinit')
    c = p.cinit;
else
    c = zeros(size(p.eb,2),1); % assume that x,y,z error coeffs are stacked end-on-end
end
    
fprintf('Starting TrACR...\n')

tti=tic; % time it

% do an initial recon of the image with the nominal kspace trajectory
G = buildG(p,c);
img = imgupdate(data,p,G);

% calculate initial cost:
cost = cost_eval(data,img,p,G);
fprintf('Starting cost: %d ... \n', cost);

%save cost and img
csv(:,:,1) = c;
imgsv(:,:,:,1) = img;
keepgoing = 1;
while keepgoing % start alternation minimization loop
    
    fprintf('Working on iter # %d ...\n',length(cost))
    
    % update k-space error coefficients
    fprintf('Updating trajectory...');
    [c,cnn] = cupdate(p,data,img,c);
    fprintf(' %d trajectory iterations taken...\n',cnn);
    
    % setup NUFFT with the new error coefficients
    G = buildG(p,c);
    
    % update image
    fprintf('Updating img...\n');
    img = imgupdate(data,p,G);
    
    % calculate cost and check exit conditions
    cost(end+1) = cost_eval(data,img,p,G);
    fprintf('Current cost: %d ... \n', cost(end));
    if length(cost) == p.maxiters+1
        fprintf('Exiting; max number of iterations reached.\n');
        keepgoing = 0;
    end
    if (cost(end-1)-cost(end) < cost(end)/100)
        fprintf('Exiting; cost is no longer decreasing significantly.\n');
        keepgoing = 0;
    end
    
    % save progress and increment iter counter
    csv(:,end+1) = c;
    imgsv(:,:,:,end+1) = img;
    
    disp('====================================================');
    
end

p.runtime = toc(tti); %in s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to build the system matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = buildG(p,c)

% reconstruct estimated total k-space trajectory
k = p.knom + reshape(p.eb*c,size(p.knom));
% normalize it
k = 2*pi*k./(p.dim./p.fov./2)./2;
% build the NUFFT (Greengard separable Gaussian kernel)
%G = gg_nufft(k,p.dim,p.Msp,p.R);
% build the gpuNUFFT

if p.useGPU
  G = gpuNUFFT(k'/2/pi, ones(size(k, 1), 1), p.osf, p.wg, p.sw, [p.dim p.dim], []);
else
  G = NUFFT(complex(k(:, 1), k(:, 2))/2/pi, ones(size(k, 1), 1, 'single'), 1, [0, 0], [p.dim p.dim], 2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image update function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img = imgupdate(data,p,G)

if strcmp(p.reconalg,'SENSE')
    
    %reconstruct images with current k-space using SENSE
    imginit = zeros([p.dim p.dim], 'single');
    img = recon_SENS(data,imginit,p.niteri,p.wi,G,p.SENSEmap);
    % replicate body coil image to all Rx coils
    img = reshape(p.SENSEmap,[p.dim,p.dim,p.ncoils]).*...
        repmat(img,[1 1 p.ncoils]);
    
else
    
    %reconstruct images with current k-space using SPIRiT
    imginit = zeros([p.dim p.dim p.ncoils], 'single');
    img = recon_SPIR(data,imginit,G,p.SPIRiTop,p.niteri,p.lam,p.wi);
    
end
% collapse the image's spatial dims
%img = permute(img,[3 1 2]);img = img(:,:).';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k-space error coefficient update function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,nn] = cupdate(p,data,img,c)

g = [];

thresh = 1e-6; % stopping threshold (in k-space traj units, eg cm^-1)

for nn = 1:p.niterk
    
    % calculate gradient
    gold = g;
    g = kgradcalc(p,data,img,c);
    
    % calculate a search direction
    if nn == 1
        dir = -g;
    else
        gamma = real(g(:)'*(g(:)-gold(:)))/real(gold(:)'*gold(:));
        dir = -g + gamma*dir;
    end
    
    % calculate a step size that reduces cost
    t = kstepcalc(p,data,img,c,g,dir,thresh);
    c = c + t*dir; % take step
    
    % exit if step is insignificant
    if t*max(abs(p.eb*dir)) < thresh ; break; end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate k-space basis derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = kgradcalc(p,data,img,c)

% p contains base k-space and error basis functions
% data contains k-space data
% img contains most recent image guess
% c contains most recent error coefficients guess

% setup NUFFT with the newest coefficients
G = buildG(p,c);

%Nc = size(data,2);
% calculate gradient
g = 0;
%if p.useparfor
%    parfor ii = 1:Nc % loop over Rx coils
%        res = data(:,ii) - G*img(:,ii); %calc residual
%        g = g + p.eb'*[real(conj(1i*2*pi*(G*(img(:,ii).*p.x(:)))).*(p.wi.*res)); ...
%            real(conj(1i*2*pi*(G*(img(:,ii).*p.y(:)))).*(p.wi.*res))];
%    end
%else
    %for ii = 1:Nc % loop over Rx coils
        res = data - G*img; %calc residual
        g = g + p.eb'*double([sum(real(conj(1i*2*pi*(G*(img.*p.x))).*(p.wi.*res)), 2); ...
            sum(real(conj(1i*2*pi*(G*(img.*p.y))).*(p.wi.*res)), 2)]);
    %end
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = kstepcalc(p,data,img,c,g,dir,thresh)

% setup NUFFT with the newest coefficients
G = buildG(p,c);

% calculate starting cost
cost = cost_eval(data,img,p,G);

% set max step size so that the maximum the traj can change in 1 update is less than
% 1/fov
tmax = min([1,1/max(p.fov)/max(abs(p.eb*dir))]);

if ~p.bls % do our own uniform backtracking line search
    
    keepgoing = 1;
    t = 0; % in case we never enter the loop
    nsteps = 10; % # of steps from 0 to max step size to consider
    
    while keepgoing && tmax*max(abs(p.eb*dir)) > thresh
        
        G = buildG(p,c + dir*tmax/nsteps);
        costt = cost_eval(data,img,p,G); %get cost w/current min step;
        
        % check whether alphamax/nsteps yields a lower cost than we currently have
        if costt < cost
            
            % if it does, go ahead and evaluate cost over current interval
            costt = [];
            %if p.useparfor
            %    parfor ii = 1:nsteps
            %        % set up NUFFT and evaluate cost with this step size
            %        G = buildG(p,c + dir*ii*tmax/nsteps);
            %        costt(ii) = cost_eval(data,img,p,G);
            %    end
            %else
                for ii = 1:nsteps
                    % set up NUFFT and evaluate cost with this step size
                    G = buildG(p,c + dir*ii*tmax/nsteps);
                    costt(ii) = cost_eval(data,img,p,G);
                end
            %end
            
            % find lowest cost and exit
            [~,mini] = min(costt);
            t = tmax/nsteps*mini;
            keepgoing = 0;
            
        else
            
            % if it doesn't, decrease alphamax and look over smaller interval
            tmax = tmax/nsteps;
            
        end
    end
    
else % use boyd's backtracking line search, which usually requires fewer cost evaluations
    
    % line search to get step
    costt = cost;
    a = 0.5; b = 0.5; t = tmax/b;
    while (costt > cost + a*t*real(g(:)'*dir)) && t*max(abs(p.eb*dir)) > thresh
        
        % reduce t
        t = b*t;
        
        % get test point
        ct = c + t*dir;
        
        % set up NUFFT and evaluate cost with this step size
        G = buildG(p,ct);
        costt = cost_eval(data,img,p,G);
        
    end
    
    if t == tmax/b || t*max(abs(p.eb*dir)) < thresh % loop was never entered or step was too small; return zero step
        t = 0;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate cost of data term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = cost_eval(data,img,p,G)

% data fidelity term (same for SPIRiT and SENSE)
%cost = 0;
%if p.useparfor
%    parfor ii = 1:p.ncoils
%        res = data(:,ii) - G*img(:,ii);
%        cost = cost + 1/2*real(res'*(p.wi.*res));
%    end
%else
    %for ii = 1:p.ncoils
        res = data - G*img;
        cost = 1/2*real(res(:)'*col(p.wi.*res));
    %end
%end

% add on the SPIRiT regularization
if strcmp(p.reconalg,'SPIRiT')
    cost = cost + p.lam*sum(col(abs(p.SPIRiTop*reshape(img,[p.dim p.dim p.ncoils]))).^2);
end
