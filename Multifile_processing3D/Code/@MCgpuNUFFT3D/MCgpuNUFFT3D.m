function  res = MCgpuNUFFT3D(k,w,b1,nz)

% Multicoil NUFFT operator for radial SOS
% using the gpuNUFFT toolbox from Andreas Schwarz
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Andreas Wetscherek, 2020
    
Nd = size(b1(:,:,1));
wg = 5; % there is no real equivalent for this parameter:
%Jd = [6,6];
osf = 1.5; % this is equivalent to the following line:
%Kd = floor([Nd*1.5]);
sw = 18; % this doesn't seem to matter too much ...
%n_shift = Nd/2;
for tt=1:size(k,3),
	kk=k(:,:,tt);
	%om = [real(kk(:)), imag(kk(:))]*2*pi;
    om = [real(kk(:)), imag(kk(:))]'; % gpuNUFFT constructor expects k-space without 2 pi and '
    % note that w can be included in gpuNUFFT, so we don't have to apply it
    % in mtimes again ...
	res.st{tt} = gpuNUFFT(om, w(:, :, tt), osf, wg, sw, Nd, [], true, true, true);
end
res.adjoint = 0;
res.imSize = size(b1(:,:,:, 1));
res.imSize2 = [size(k,1),size(k,1),nz];
res.dataSize = size(k);
res.w = sqrt(w);
res.b1 = b1;

res = class(res,'MCgpuNUFFT3D');

