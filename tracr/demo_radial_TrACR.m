% Demo Trajectory Auto-Compensation Reconstruction (TrACR)for radial k-space trajectories

%addpath ../util/
%addpath(genpath('../../gpuNUFFT/gpuNUFFT'));

% using SENSE image reconstructions
% load in data-------------------------------------------------------------
%fprintf('Loading data...\n');
%load('test_radial_brain_multicoil.mat');
npts = size(knom,1); % number of points per radial line
nproj = size(knom,2); % number of radial projections

% TrACR parameters---------------------------------------------------------
p.data = single(data);   % k-space data in: [(# readout pts * shots) X # coils]
p.ncoils = size(data,2); % # Rx coils
p.knom = single([real(knom(:)) imag(knom(:))]);   % nominal k-space trajectory form, cycles/cm [# points per projection x # projections]
p.usedc = 1;     % use density compensation
p.dim = par.kx_range(2) + 1;     % matrix size, 1D (2nd assumed equal)
p.fov = par.FOVx / 10;    % field of view in cm
p.maxiters = 100;  % maximum # of outer iterations to run TrACR
p.niteri = 2;    % # of iterations per image update
p.niterk = 5;    % # of iterations per kspace update
p.imgcart = single(imgcart);    % Cartesian image from which to derive coil sensitivities
p.reconalg = 'SENSE'; % use SENSE image reconstruction (TrACR can also do SPIRiT if the correct parameters are passed)

%NUFFT parameters----------------------------------------------------------
%p.Msp = 6; % spreading parameter
%p.R = 2; % oversampling ratio

if p.useGPU
  %gpuNUFFT parameters
  p.osf = 2; % oversampling: 2 1.5 1.25
  p.wg = 7; % kernel width: 3 5 7
  p.sw = 8; % parallel sectors' width: 8 12 16
end

% set up image grid locations - used in k-space derivative calculation
[p.x,p.y] = ndgrid(-p.fov/2:p.fov/p.dim:p.fov/2-p.fov/p.dim);
%[x,y] = ndgrid(-p.fov/2:p.fov/p.dim:p.fov/2-p.fov/p.dim);
%p.x = x(:);
%p.y = y(:);

% density compensation weights: for GA radial, approximate as |kx+1i*ky|
if p.usedc
    %p.wi = abs(p.knom(:,1) + 1i*p.knom(:,2))./max(abs(p.knom(:,1) + 1i*p.knom(:,2)));
    %p.wi(isnan(p.wi))=1;
    %p.wi(isinf(p.wi))=1;
    p.wi = w;
else
    p.wi = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up k-space error basis (for radial)
% Assumption: The k-space errors comprise shifts of each line in kx-ky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Setting up error basis...\n');
eb = kron(speye(nproj),ones(npts,1));
p.eb = [eb 0*eb;0*eb eb]; % k-traj error basis (npts*nproj*(kx + ky) x nproj)

%%%%%%%%%%%%%%%%%%%%%%%%
% estimate sensitivites
%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Estimating Sensitivities...\n');
sosimg = sqrt(sum(abs(p.imgcart.^2),3));
%sens = reshape(double(p.imgcart./repmat(sosimg,[1 1 p.ncoils])), [p.dim*p.dim, p.ncoils]);
%p.SENSEmap = double(sens);
p.SENSEmap = reshape(p.imgcart./repmat(sosimg,[1 1 p.ncoils]), [p.dim*p.dim, p.ncoils]);

% run TrACR----------------------------------------------------------------
%parfor ii=1:feature('numCores')
%    nufftdummy = gg_nufft(zeros([2 2]),3,6,2);
%end
[img,c,cost,imgsv,csv,p] = TrACR(p.data,p);
k = p.knom + reshape(p.eb*c,size(p.knom)); % corrected trajectory (first col: kx, second col: ky)

fprintf(['TrACR finished with: %4.0f iterations. \n ', ...
    'Total compute time: %3.0f minutes %3.0f seconds\n'], ...
    size(cost,2)-1,p.runtime/60, rem(p.runtime,60));

%display results-----------------------------------------------------------
%img_initial = imgsv(:,:,:,1);
%img_initial = sqrt(sum(abs(img_initial).^2,3)).';
%img_final = reshape(img,p.dim,p.dim,[]);
%img_final = sqrt(sum(abs(img_final).^2,3)).';

%figure; subplot(2,2,1); imagesc(img_initial,[0 max([img_initial(:); img_final(:)])]);
%colormap gray; title('initial sum-of-squares image'); axis off; axis image;
%subplot(2,2,2); imagesc(img_final,[0 max([img_initial(:); img_final(:)])]);
%colormap gray; title('final sum-of-squares image'); axis off; axis image;
%subplot(2,2,3); imagesc(img_final-img_initial);
%colormap gray; title('difference'); axis off; axis image;
%subplot(2,2,4); semilogy(cost); xlabel('iteration #'); ylabel('data cost term');
%title('cost');

