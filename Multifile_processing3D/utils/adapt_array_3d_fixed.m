function [recon,cmap]=adapt_array_3d(yn,rn,norm);

% Reconstruction of array data and computation of coil sensitivities based 
% on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
%-------------------------------------------------------------------------
%	Input:
%	yn: array data to be combined [ny, nx, nz, nc]. 
%	rn: data covariance matrix [nc, nc].
%	norm: =1, normalize image intensity
%
%	Output:
%	recon: reconstructed image [ny, nx, nz].
%	cmap: estimated coil sensitivity maps [ny, nx, nz, nc].
%--------------------------------------------------------------------------
% Ricardo Otazo
% CBI, New York University
%--------------------------------------------------------------------------
% made 3D by Andreas Wetscherek

yn=permute(yn,[4,1,2,3]);
[nc,ny,nx,nz]=size(yn);
if nargin<3, norm=0; end
if nargin<2, rn=eye(nc);end

% find coil with maximum intensity for phase correction
[mm,maxcoil]=max(sum(sum(sum(permute(abs(yn),[4 3 2 1])))));

bs1=4;  %x-block size
bs2=4;  %y-block size
bs3=2;  %z-block size
st=1;   %increase to set interpolation step size (do we have to ?)
% this code is anyway not optimized for speed :P

% not sure why only ny had the round() around it ...
wsmall=zeros(nc,round(ny./st),round(nx./st), round(nz./st));
cmapsmall=zeros(nc,round(ny./st),round(nx./st), round(nz./st));

for z=st:st:nz
for x=st:st:nx
for y=st:st:ny
    %Collect block for calculation of blockwise values
    ymin1=max([y-bs1./2 1]);                   
    xmin1=max([x-bs2./2 1]);  
    zmin1=max([z-bs3./2 1]);
    % Cropping edges
    ymax1=min([y+bs1./2 ny]);                 
    xmax1=min([x+bs2./2 nx]); 
    zmax1=min([z+bs3./2 nz]);              

    ly1=length(ymin1:ymax1);
    lx1=length(xmin1:xmax1);
    lz1=length(zmin1:zmax1);
    m1=reshape(yn(:,ymin1:ymax1,xmin1:xmax1,zmin1:zmax1),nc,lx1*ly1*lz1);
      
    m=m1*m1'; %signal covariance
      
    % eignevector with max eigenvalue for optimal combination
    [e,v]=eig(inv(rn)*m);                    
                                               
    v=diag(v);
    [mv,ind]=max(v);
      
    mf=e(:,ind);                      
    mf=mf/(mf'*inv(rn)*mf);               
    normmf=e(:,ind);
    
    % Phase correction based on coil with max intensity
    mf=mf.*exp(-j*angle(mf(maxcoil)));        
    normmf=normmf.*exp(-j*angle(normmf(maxcoil)));

    wsmall(:,y./st,x./st,z./st)=mf;
    cmapsmall(:,y./st,x./st,z./st)=normmf;
end
end
end

recon=zeros(ny,nx,nz);

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude 
% pixels between +1 and -1 pixels.
wfull = conj(wsmall);
cmap = cmapsmall;

recon=squeeze(sum(wfull.*yn));   %Combine coil signals. 
% normalization proposed in the abstract by Griswold et al.
if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

cmap=permute(cmap,[2,3,4,1]);