function [img,flag,relres,iter,resvec,lsvec] = recon_SENS(data,imginit,niteri,wi,G,SENSEmap,tol)
% reconstruct an image from non-Cartesian acquisitions, using SENSE

%data       =   initial k-space data (must be a double!)
%imginit    =   initial image
%niteri     =   number of lsqr iterations
%wi         =   density weights
%G          =   NUFFT object (set with Gnufft)
%SENSEmap   =   map of coil sensitivities (must be a double!)

%SENSEmap = SENSEmap(:);
N = numel(imginit);     %number of image points
imSize = size(imginit); %image size
dataSize = size(data);  %k-space data size

wi=sqrt(wi(:));
b = col(data.*wi);
if ~exist('tol','var')
    tol=1e-6;
end
%perform least-squares recon
[img,flag,relres,iter,resvec,lsvec] = lsqr(@(x,transp_flag)afun(x,wi,dataSize,imSize,G,SENSEmap,transp_flag), b, tol, niteri,[],[], imginit(:));

img = reshape(img,imSize);

function [y, transp_flag] = afun(x,wi,dataSize,imSize,G,SENSEmap,transp_flag)

if strcmp(transp_flag,'transp')

    %reshape to orig data size
    %x1 = reshape(x(1:prod(dataSize)),dataSize);
    x1 = reshape(x,dataSize);
    %numcoils = dataSize(end);
    %y = zeros([imSize, numcoils]);
    %parfor ii = 1:numcoils
    %for ii = 1:numcoils
        %type I NUFFT
        y = reshape(G'*(wi.*x1), size(SENSEmap));
    %end
    %y = y(:);
    y = conj(SENSEmap).*y;
    %y = reshape(y,[prod(imSize),dataSize(end)]);
    y = sum(y,2);
    %y = double(y);

else

    xt=reshape(SENSEmap.*x, [imSize, dataSize(end)]);
    %xt=reshape(xt,[imSize dataSize(end)]);
    %numcoils=dataSize(end);
    %y1=zeros(dataSize);
    %parfor ii = 1:numcoils
    %for ii = 1:numcoils
        %xtall = xt(:,:,ii);
        %type II NUFFT
        y=wi.*(G*xt);
        y=y(:);
    %end
    %y = y1(:);

end
