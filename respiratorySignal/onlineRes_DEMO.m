%%
clear; clc
addpath('../../../arrShow/');
addpath('./Dataset');
addpath('./util')

%parameters
initial = 150;  %how many initial spokes we use for getting the eigenvectors
thresh = 0.95;  %threshold for coilClustering function
online_plot = 1;%1->display online plot 
                %0->don't display the plot (FASTER)

load kc.mat                %contains k-space center (kc) data with the dims: [nkz nSpokes ncoils]
[nz,nSpokes_full,nc]=size(kc);  %with the dims: [nkz nSpokes ncoils]
nspokes_online = size(kc,2) - initial;

if online_plot == 1
    %getting the reference/correct res. signal with the original XD-GRASP code
    Res_Signal_correct = XDGRASP_res_sig(kc,thresh);

    %online plot init:
    plot(-Res_Signal_correct);      %background plot for reference
    hold on
    Res_Signal = zeros(nSpokes_full,1);
    xs = 1:length(Res_Signal);
    p = plot(xs,Res_Signal);
    p.YDataSource = 'Res_Signal';
end

% First we get the initial data and construct the PCA eigenvectors
    %we want to get normalization values, vectors PC1 and PC2 for each coil
    %Res_sig and cluster (from CoilClustering)

kc_init = kc(:,1:initial,:);
[nz,nSpokes,nc]=size(kc_init);
ZIP=abs(ifft(kc_init,400,1));

%Normalization of each projection in each coil element (I don't think that this is necessary, the code should work without)
for ii=1:nc
    for jj=1:nSpokes
        max_coil_spoke=max(ZIP(:,jj,ii));
        min_coil_spoke=min(ZIP(:,jj,ii));
        ZIP(:,jj,ii)=(ZIP(:,jj,ii)-min_coil_spoke)./(max_coil_spoke-min_coil_spoke);
    end
end


 kk = 1; clear PCs   
 V1_all = zeros(400, nc);   %eigenvectors with the 1st highest ev for all coils
 V2_all = zeros(400, nc);   %eigenvectors with the 2nd highest ev for all coils
 MIN = zeros(length(nc),2); %for saving min values for online normalization
 MAX = zeros(length(nc),2); %for saving maxs values for online normalization
 
 for ii=1:nc    %loop through the coils
    ZIP_coil = squeeze(ZIP(:,:,ii))'; %kspace data after fft in kz, only for one coil
    covariance = cov(ZIP_coil); %gives the covariant matrix, diagonal entries are variance
                                %var in the first dim, nondiagonal entries are either correlation or
                                %or anti correlation
    [V E]= eig(covariance);     %V eigenvectors and E eigenvalues of my covariant matrix 
                                %!eigenvector with the largest eigenvalue is the direction of the
                                %maximum variance of the data -> PC1
    eigenvalues = diag(E);
    [sorted order] = sort(eigenvalues, 'descend');
    eigenvalues = sorted;
    
    V_ordered = V(:,order); %order the eigenvectors in the eigenvalues order from the vector with highest to lowest eigv:
    V1_all(:,ii) = V_ordered(:,1); %save the PC1 vector for each coil
    V2_all(:,ii) = V_ordered(:,2); %save the PC2 vector for each coil
    
    %projection (dot prod) of the 2D kcdata onto the principal vectors
    PC = (V_ordered' * ZIP_coil')'; % for every spoke for every coil
                                    % we have projected our data onto vectors with most variance

    % Take the first two principal components from each coil element.
     for jj=1:2
         %we throw away all the other components besides the projections
         %on the vector PC1 and PC2 because on those vectors we should see 
         %the biggest variance of the data
         PC_smoothed = smooth(PC(:,jj),6, 'lowess'); %moving average smoothing
        
         %normalize
         MIN(ii,jj) = min(PC_smoothed);             % we are saving the max and min values for each principal component for the future online normalization 
         PC_smoothed_shifted = PC_smoothed - MIN(ii,jj); 
         MAX(ii,jj) = max(PC_smoothed_shifted);     % we are saving the max and min values for each principal component for the future online normalization
         PC_normalized = PC_smoothed_shifted ./ MAX(ii,jj) ;

         PCs(:,kk) = PC_normalized;
         kk = kk + 1;
     end    
 end

%coil clustering to find the res. motion signal
[Res_Signal_init, global_clustering] = CoilClustering(PCs, thresh);
    
%Normalize the res. signal
min_res = min(Res_Signal_init(:));
Res_Signal_init=Res_Signal_init-min_res;    %moving the minimum to zero
max_res = max(Res_Signal_init(:));
Res_Signal_init=Res_Signal_init./max_res;   %moving the max to one 
Res_Signal(1:length(Res_Signal_init)) = Res_Signal_init;    
%get global clustering for PC1 and PC2 separately (we are using this for all the other future spokes -> 'global'):
clusteringPC1 = global_clustering(1:2:length(global_clustering)); %either add or substract the signal from the coil or don't use it at all
clusteringPC2 = global_clustering(2:2:length(global_clustering)); %either add or substract the signal from the coil or don't use it at all
coils_combined = sum(abs(global_clustering));   %how many coils are combined -> for normalization  


tic; %timing
%Now perform the same projection with the other spokes - 'online res. signal acquisition'
for spoke = 151:size(kc,2)
    kc_spoke = squeeze(kc(:,spoke,:));
    ZIP_spoke=abs(ifft(kc_spoke,400,1));
    
    %normalization within the spoke
    for ii=1:nc
            max_coil_spoke=max(ZIP_spoke(:,ii));
            min_coil_spoke=min(ZIP_spoke(:,ii));
            ZIP_spoke(:,ii)=(ZIP_spoke(:,ii)-min_coil_spoke)./(max_coil_spoke-min_coil_spoke);
    end

    
    SIG = 0;        %the scalar value of res. signal from each spoke (normalized)
    for ii = 1:nc   %loop over the coils               ( this is possible to rewrite without a loop with mtimesx package  mtimesx(V1_all(:,ii)',ZIP_spoke) )
        SIG1_iter = V1_all(:,ii)'*ZIP_spoke(:,ii);      %projection onto the eigenvector PC1
        SIG1_iter = (SIG1_iter - MIN(ii,1))./MAX(ii,1); %normalization?
        SIG2_iter = V2_all(:,ii)'*ZIP_spoke(:,ii);      %projection onto the eigenvector PC2
        SIG2_iter = (SIG2_iter - MIN(ii,2))./MAX(ii,2); %normalization? 
        SIG = SIG + clusteringPC1(ii)*SIG1_iter + clusteringPC2(ii)*SIG2_iter; %we sum all the contributions
    end

    SIG = SIG/coils_combined;       %normalization of the online 'coil clustering'
    SIG = (SIG-min_res)/max_res;    %final normalization
    Res_Signal(spoke) = SIG;
    
    %plot the resulting (updated) signal online each time iteration:
    if online_plot ==1
        refreshdata
        drawnow
        %pause(0.1)
    end
end    

%timing info
t = toc; 
disp(['avg. time per spoke: ',num2str(t/nspokes_online)]); %it makes difference if plotting online or without online plotting
%show the resulting plot
    %plot(Res_Signal);
    %title('Online acquired res. signal');