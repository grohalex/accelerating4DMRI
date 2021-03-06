%% res signal for changing number of spokes
% This program finds the angles between the eigenvector reconstructed from
% k-space center (kc), containing all the spokes AND eigenvectors
% reconstructed from smaller subsets (containing smaller number of spokes) 
% of kc. It compares eigenvectors from PCA (Principal Component Analysis),
% PC1 and PC2  
%   -> the conlusion is that the vectors are well aligned even
%   for smaller number of spokes but you can see that they can be flipped by 
%   pi radians resulting in flipped res. signal which is a possible issue
%
%The code was originally developed for XD-GRASP liver MRI 
%as described in:
%  
%  Feng L, Axel L, Chandarana H, Block KT, Sodickson DK, Otazo R. 
%  XD-GRASP: Golden-angle radial MRI with reconstruction of 
%  extra motion-state dimensions using compressed sensing
%  Magn Reson Med. 2016 Feb;75(2):775-88
% (c) Li Feng, 2016, New York University
%
% The code was modified by Alexander Groh, 2020, ICR

clear; clc
addpath('../../arrShow/');
addpath('./Dataset');
addpath('./Code');
%loading the center of kspace data (kx = ky = 0) for each kz slice
load kc.mat; %load kc125; %load kc115 

compare_spokes = flip([10:50]); %[713 650 600 550 500 450 400 355 142 100 71 50 35 ]; %spokes for different recons which will be compared
Res_Signal = zeros(length(compare_spokes), compare_spokes(1));
VectorsPC1 = zeros(400,length(compare_spokes));
VectorsPC2 = zeros(400,length(compare_spokes));
for iter = 1:length(compare_spokes)
    clearvars -except iter compare_spokes kc Res_Signal VectorsPC1 VectorsPC2
    
    kc_iter = kc(:, 1:compare_spokes(iter), :); %we define only a subset of kc (only limited number of spokes)
    [nz,nSpokes,nc]=size(kc_iter);
    
    %nz: number of slice
    %ntviews: number of acquired spokes
    %nc: number of coil elements

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Respiratory motion detection

    % Generate the z-projection profiles
    %ZIP: Projection profiles along the Z dimension with interpolation.
    ZIP=abs(ifft(kc_iter,400,1));

    %Normalization of each projection in each coil element
    for ii=1:nc
        for jj=1:nSpokes
            max_coil_spoke=max(ZIP(:,jj,ii));
            min_coil_spoke=min(ZIP(:,jj,ii));
            ZIP(:,jj,ii)=(ZIP(:,jj,ii)-min_coil_spoke)./(max_coil_spoke-min_coil_spoke);
        end
    end

    %as(ZIP)
    %% Perform PCA on each coil element (principal component analysis)

    kk=1;clear PCs
    for ii=1:nc % loop through the coils
        ZIP_permuted = permute(ZIP(:,:, ii),[1 3 2]); 
        ZIP_reshaped = (squeeze(ZIP_permuted))'; %dims: [fft_kz*1, spokes(time) ]'
      
        covariance = cov(ZIP_reshaped); %gives the covariant matrix, diagonal entries are variance
        %var in the first dim, nondiagonal entries are either correlation or
        % or anti correlation

        [V E]= eig(covariance); % I get V eigenvectors and E eigenvalues of my covariant matrix 
        %!eigenvector with the largest eigenvalue is the direction of the
        %maximum variance of the data -> PC1!!!

        eigenvalues = diag(E);
        [sorted order] = sort(eigenvalues, 'descend');
        eigenvalues = sorted;
        %order the eigenvectors in the eigenvalues order from the vector with highest to lowest eigv:
        V_ordered = V(:,order); 
        %projection (dot prod) of the 2D kcdata onto the principal vectors
        PC = (V_ordered' * ZIP_reshaped')'; % for every spoke for every coil
        %so now we have projected our data onto vectors with most variance

        % Take the first two principal components from each coil element.
         for ii=1:2
             %we throw away all the other components besides the projections
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

    %coil clustering to find the res. motion signal
    thresh = 0.95;
    [Res_Signal_iter, cluster] = CoilClustering(PCs, thresh);%CoilClustering(PCs, thresh);

    %Normalize the res. signal
    Res_Signal_iter=Res_Signal_iter-min(Res_Signal_iter(:)); %moving the minimum zero
    Res_Signal_iter=Res_Signal_iter./max(Res_Signal_iter(:)); %moving the max to one 

    Res_Signal( iter, 1:length(Res_Signal_iter)) = Res_Signal_iter; 
    % Plot the respiratory motion signal on the projection profiles
    %{
    imagesc(abs(ZIP(:,:,5))),axis image,colormap(gray), axis off,title('Respiratory Motion from rewritten XD-GRASP')
    hold on
    plot(-Res_Signal(:)*100+220,'r')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    hold off
    %}
    VectorsPC1(:,iter) = V_ordered(:,1);
    VectorsPC2(:,iter) = V_ordered(:,2);
end    
%find the angles between the original PC1 vector and PC1 vectors form
%smaller kc
for i = 1:length(compare_spokes) %find the angles between the vectors PC1:  cos(a) = V.W/(||V||*||W||)
    angles(i) = acos(dot(VectorsPC1(:,1),VectorsPC1(:,i))/(norm(VectorsPC1(:,1))*norm(VectorsPC1(:,i))));
end

%find the angles between the original PC2 vector and PC2 vectors from
%smaller kc
for i = 1:length(compare_spokes) %find the angles between the vectors PC2:  cos(a) = V.W/(||V||*||W||)
    angles2(i) = acos(dot(VectorsPC2(:,1),VectorsPC2(:,i))/(norm(VectorsPC2(:,1))*norm(VectorsPC2(:,i))));
end
%plot all the angles (projections)
plot(1:length(compare_spokes), angles)
 title('Changing angle between the eigenvectors obtained with different numbers of spokes')
 xlabel(num2str([compare_spokes]))
 hold on 
 plot(angles2)
 legend('PC1 eigenvectors','PC2 eigenvectors')