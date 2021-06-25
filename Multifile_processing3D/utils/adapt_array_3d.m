function [recon,cmap,wfull]=adapt_array_3d(yn,rn,norm);

% Reconstruction of array data and computation of coil sensitivities based 
% on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
%-------------------------------------------------------------------------
%	Input:
%	yn: array data to be combined [ny, nx, nc]. 
%	rn: data covariance matrix [nc, nc].
%	norm: =1, normalize image intensity
%
%	Output:
%	recon: reconstructed image [ny, nx].
%	cmap: estimated coil sensitivity maps [ny, nx, nc].
%--------------------------------------------------------------------------
% Ricardo Otazo
% CBI, New York University
%--------------------------------------------------------------------------
%
% modified by Andreas Wetscherek for 3D (a.wetscherek@icr.ac.uk) 07/07/2020
yn = permute(yn, [4, 1, 2, 3]);
[nc,ny,nx,nz]=size(yn);
if nargin<3, norm=0; end
if nargin<2, rn=eye(nc);end

rni = inv(rn); %#ok<*MINV>

% find coil with maximum intensity for phase correction
[~,maxcoil]=max(sum(sum(sum(abs(yn), 4), 3), 2));   

bs=[8,8,4];  %x, y, z block size
st=bs/2;     %x, y, z interpolation step size

nx_low = round((nx-bs(1)+st(1))/st(1));
ny_low = round((ny-bs(2)+st(2))/st(2));
nz_low = round((nz-bs(3)+st(3))/st(3));

wsmall=zeros(nc, nx_low+2, ny_low+2, nz_low+2);
cmapsmall=zeros(nc, nx_low+2, ny_low+2, nz_low+2);

for x=st(1):st(1):nx-1
for y=st(2):st(2):ny-1
for z=st(3):st(3):nz-st(3)
    %Collect block for calculation of blockwise values                   
    xmin=x-st(1)+1;
    ymin=y-st(2)+1;
    zmin=z-st(3)+1;
    
    xmax=x-st(1)+bs(1);   
    ymax=y-st(2)+bs(2);               
    zmax=z-st(3)+bs(3);                 

    m1=reshape(yn(:,ymin:ymax,xmin:xmax,zmin:zmax),nc,prod(bs));
      
    m=m1*m1'; %signal covariance
      
    % eignevector with max eigenvalue for optimal combination
    [e,v]=eig(rni*m);                     
                                               
    v=diag(v);
    [~,ind]=max(v);
      
    mf=e(:,ind);                      
    mf=mf/(mf'*rni*mf);               
    normmf=e(:,ind);
    
    % Phase correction based on coil with max intensity
    mf=mf.*exp(-1i*angle(mf(maxcoil)));        
    normmf=normmf.*exp(-1i*angle(normmf(maxcoil)));

    wsmall(:,y/st(2)+1,x/st(1)+1,z/st(3)+1)=mf;
    cmapsmall(:,y/st(2)+1,x/st(1)+1,z/st(3)+1)=normmf;
end
    wsmall(:, :, :, 1) = wsmall(:, :, :, 2);
    wsmall(:, :, :, end) = wsmall(:, :, :, end-1);
    cmapsmall(:, :, :, 1) = cmapsmall(:, :, :, 2);
    cmapsmall(:, :, :, end) = cmapsmall(:, :, :, end-1);
end
    wsmall(:, 1, :, :) = wsmall(:, 2, :, :);
    wsmall(:, end, :, :) = wsmall(:, end-1, :, :);
    cmapsmall(:, 1, :, :) = cmapsmall(:, 2, :, :);
    cmapsmall(:, end, :, :) = cmapsmall(:, end-1, :, :);
end

wsmall(:, :, 1, :) = wsmall(:, :, 2, :);
wsmall(:, :, end, :) = wsmall(:, :, end-1, :);
cmapsmall(:, :, 1, :) = cmapsmall(:, :, 2, :);
cmapsmall(:, :, end, :) = cmapsmall(:, :, end-1, :);

[Xq, Yq, Zq] = meshgrid(((1:nx) - (bs(1) + 1)/2) / st(1) + 2, ...
    ((1:ny) - (bs(2) + 1)/2) / st(2) + 2, ((1:nz) - (bs(3) + 1)/2) / st(3) + 2);

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude 
% pixels between +1 and -1 pixels.
for i=1:nc
    wfull(i,:,:,:)=interp3(squeeze(abs(wsmall(i,:,:,:))),Xq,Yq,Zq,'linear').*exp(1i.*angle(interp3(squeeze(wsmall(i,:,:,:)),Xq,Yq,Zq,'linear')));
    cmap(:,:,:,i)=interp3(squeeze(abs(cmapsmall(i,:,:, :))),Xq,Yq,Zq,'linear').*exp(1i.*angle(interp3(squeeze(cmapsmall(i,:,:, :)),Xq,Yq,Zq,'linear')));
end
recon=squeeze(sum(conj(wfull).*yn));   %Combine coil signals. 
% normalization proposed in the abstract by Griswold et al.
if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

wfull = permute(wfull, [2, 3, 4, 1]);
