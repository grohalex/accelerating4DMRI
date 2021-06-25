function  res = MCNUFFT_GPU3D(k_u,w_u,b1,nz)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng, 2016
%
% modified by Andreas Wetscherek, 2020, The Institute of Cancer Research
% to perform computations on GPU

% The basic idea is to just perform a couple of matrix multiplications and
% keep all data on the GPU at all times. The gpuNUFFT MATLAB interface does
% not allow for this, as data needs to be in normal RAM for wrapping. The
% speedup gained (hopefully) from the implementation at hand comes from not
% having any transfer costs.
% However, the complete (dense) Matrix of Fourier coefficients is too large
% to store on GPU (at least my 48 GB are not enough). But it is possible to
% break it up along the dimension of image positions, as we're on a
% Cartesian grid. If we move by a vector (xoffs, yoffs), we basically can
% reuse the Fourier coefficient matrix for (0, 0) and just insert another
% multiplication with a diagonal matrix given by:
% exp(2pi * 1i * om \cdot (xoffs, yoffs))

% we want to make this (dense) Matrix of Fourier coefficients as large as
% possible, so the limitation here is the number of elements, but we need a
% separate Matrix for each respiratory phase. So we might run into memory
% limitations after all...

gpu = gpuDevice;
memavail = gpu.AvailableMemory * 0.9;

Nd = size(b1(:,:,1));
%Jd = [6,6];
%Kd = floor([Nd*1.5]);
n_shift = Nd/2;
[x, y] = ndgrid(-n_shift(1):Nd(1)-n_shift(1)-1);

% find a useful number to divide the (dense) Matrix of fourier coefficients
% along each dimension, such that the resulting smaller matrices for all
% respiratory phases still fit in memory
n_submatrix = single(ceil(numel(x) * size(k_u, 1) * size(k_u, 2) / intmax('int32') / 2) * 2);
while mod(size(x, 2), n_submatrix*n_submatrix) ~= 0 ...
        || numel(x) * size(k_u, 1) * size(k_u, 2) * size(k_u,3) * 8 / n_submatrix > memavail
    n_submatrix = n_submatrix + 1;
end
n_posperdim = size(x, 1) / n_submatrix;

for tt=1:size(k_u,3),
	kk=k_u(:,:,tt);
	om = [real(kk(:)), imag(kk(:))]*2*pi;
	%res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
    res.st{tt} = exp(1i * (om(:, 1).' .* col(x(1:n_posperdim,1:n_posperdim) + n_shift(1)) ...
                         + om(:, 2).' .* col(y(1:n_posperdim,1:n_posperdim) + n_shift(2))));
    for ss = 1:n_submatrix
        for rr = 1:n_submatrix
            res.sp{tt,ss,rr} = exp(1i * (om(:, 1) .* x(n_posperdim*(ss-1)+1,1) ...
                                       + om(:, 2) .* y(1,n_posperdim*(rr-1)+1)));
        end
    end
end

res.adjoint = 0;
res.imSize = size(b1(:,:,:,1));
% this did not make really sense to me... but keeping it for numerical
% comparison. The second value should be size(k_u,2):
res.imSize2 = [size(k_u,1),size(k_u,1),nz]; 
res.dataSize = size(k_u);
res.w = sqrt(w_u);
res.b1 = b1;
res = class(res,'MCNUFFT_GPU3D');

