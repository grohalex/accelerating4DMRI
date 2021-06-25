function  res = MCNUFFT3D(k,w,b1,pf_mask)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
% pf_mask: mask to account for sampled k-space locations
%
% Li Feng, 2016
    
Nd = size(b1(:,:,1)); % even though b1 is 4D here, this will be fine
Jd = [6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;
for tt=1:size(k,3),
	kk=k(:,:,tt);
	om = [real(kk(:)), imag(kk(:))]*2*pi;
	res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
end
res.adjoint = 0;
res.imSize = size(b1(:,:,1));
res.imSize2 = [size(k,1),size(k,1)];
res.dataSize = size(k);
res.w = sqrt(w);
res.b1 = b1;
res.pf_mask = pf_mask;
res = class(res,'MCNUFFT3D');
