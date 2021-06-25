function  res = NUFFT_irt(k,w,imSize)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng, 2016
    
Nd = imSize;
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
%for tt=1:size(k,3),
%	kk=k(:,:,tt);
	om = [real(k(:)), imag(k(:))]*2*pi;
	res.st = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
%end
res.adjoint = 0;
res.imSize = imSize;
res.imSize2 = [size(k,1),size(k,2)];
res.dataSize = size(k);
res.w = sqrt(w);
%res.b1 = b1;
res = class(res,'NUFFT_irt');

