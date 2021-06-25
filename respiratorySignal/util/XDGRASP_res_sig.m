function signal = XDGRASP_res_sig(kc, thresh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  In the following code I am using modified version of Matlab 
%  source code for XD-GRASP liver MRI  as described in:
% 
%  
%  Feng L, Axel L, Chandarana H, Block KT, Sodickson DK, Otazo R. 
%  XD-GRASP: Golden-angle radial MRI with reconstruction of 
%  extra motion-state dimensions using compressed sensing
%  Magn Reson Med. 2016 Feb;75(2):775-88
% 
%
%  The source code uses the following external packages:
%    - NUFFT toolkit by Jeffrey Fessler 
%      (http://www.eecs.umich.edu/~fessler/)
%    - Non-linear conjugate gradient algorithm by Miki Lustig
%      (http://www.eecs.berkeley.edu/~mlustig/Software.html)
%    - Coil clustering code by Tao Zhang 
%      (http://web.stanford.edu/~tzhang08/software.html)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nz,nSpokes,nc]=size(kc);
ZIP=abs(ifft(kc,400,1)); %ZIP: Projection profiles along the Z dimension with interpolation.

%Normalization of each projection in each coil element
for ii=1:nc
    for jj=1:nSpokes
        max_coil_spoke=max(ZIP(:,jj,ii));
        min_coil_spoke=min(ZIP(:,jj,ii));
        ZIP(:,jj,ii)=(ZIP(:,jj,ii)-min_coil_spoke)./(max_coil_spoke-min_coil_spoke);
    end
end

kk=1;clear PCs
for ii=1:nc % loop through the coils
    ZIP_permuted = permute(ZIP(:,:, ii),[1 3 2]);   %dims: [fft_kz 1 spokes(time)]
    ZIP_reshaped = (squeeze(ZIP_permuted))';        %dims: [fft_kz, spokes(time)]'
    %PCA
    covariance = cov(ZIP_reshaped); %gives the covariant matrix, diagonal entries are variance
    [V E]= eig(covariance);         % I get V eigenvectors and E eigenvalues of my covariant matrix 
    eigenvalues = diag(E);
    [sorted order] = sort(eigenvalues, 'descend');
    eigenvalues = sorted;
    V_ordered = V(:,order);         %order the eigenvectors in the eigenvalues order from the vector with highest to lowest eigv:
    
    %Projection (dot prod) of the 2D kcdata onto the principal vectors
    PC = (V_ordered' * ZIP_reshaped')';
    
    %Take the first two principal components from each coil element.
     for ii=1:2
         %throw away all the other components besides the projections
         %on the vector PC1 and PC2 because on those vectors we should see 
         %the biggest variance of the data
         PC_smoothed = smooth(PC(:,ii),6, 'lowess'); %moving average smoothing
         
         %normalize
         PC_smoothed_shifted = PC_smoothed - min(PC_smoothed); 
         PC_normalized = PC_smoothed_shifted ./ max(PC_smoothed_shifted);
         
         PCs(:,kk) = PC_normalized;
         kk = kk + 1;
     end    
end
%coil clustering to find the res. motion signal -> choosing/combining coils
[signal, cluster] = CoilClustering(PCs, thresh);

%Normalize the res. signal
signal=signal-min(signal(:)); %moving the minimum to zero
signal=signal./max(signal(:)); %moving the max to one 
end
