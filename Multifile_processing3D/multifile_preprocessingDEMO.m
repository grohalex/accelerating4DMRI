%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Matlab source code for XD-GRASP liver MRI 
%  as described in:
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

%  If you use this code, please cite the above publication.
%
%  (c) Li Feng, 2016, New York University
%  Li.Feng@nyumc.org
%
% Modified for 3D by Andreas Wetscherek (a.wetscherek@icr.ac.uk) 09/07/2020
% Modified by Alexander Groh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Description
% multi file preprocessing - to process one or more files and create a
% folder with a identifying name
% it also creates a table which contains all the information about every
% file that was preprocessed

clear
clc
addpath('F:\ICR_Studentship\gpuNUFFT\gpuNUFFT')%('E:/MATLAB/gpuNUFFT/gpuNUFFT'); % if gpuNUFFT is to be used
addpath(genpath('Code'));
addpath('utils');

par.useFullGPU = false; % this is not fully implemented yet ...
par.noGPU = false;
XD = false;%true; % perform XD-GRASP algorith with the preprocessed data?
table = [];
%write the path to each file that you want to preprocess in this cell (multifile - this can preprocess as many files as are put in the cell)
filenames ={'F:\ICR_Studentship\data\20191213_115503_P131\raw_000.data',
            'F:\ICR_Studentship\data\20191120_125927_P131\raw_000.data' };
%write the dicom directories, for FOV and positioning / referencefile
dicoms = {  'F:\ICR_Studentship\data\20191213_115503_P131\scans\205_DelRec - T1 3DVane w_o fatsat thorax 1_55 mm (5_38 duration)\DICOM\*.dcm',
            'F:\ICR_Studentship\data\20191120_125927_P131\scans\203_DelRec - T1 3DVane w_o fatsat thorax 1_55 mm (5_38 duration)\DICOM\*.dcm'};

N = strcat('multifile_processing_', date);
mkdir(N); cd(N);
homedir = pwd;

for File=1:length(filenames) 
    % raw data set (.list + .data) as from MR-Linac
    disp(['Preprocessing file ', num2str(File), ': ',char(filenames(File))]) %%%
    par.filename = char(filenames(File));%F:\ICR_Studentship\data\20191213_115503_P131\raw_000.data';
    a = strsplit(par.filename, '\');
    folder = char(a(numel(a)-1));

    % working directory in which results are stored that can be achieved before
    % the XDGRASP algorithm is started, e.g. density compensation, trajectory
    % correction, coil sensitivity maps ...
    par.workspace = strcat(homedir, '\preprocessing_',folder) % this created folder multiprocessing_20191213_115503_P131
    cd(homedir), addpath(homedir)
    mkdir(par.workspace);

    % matching Dicom images from MR-Linac (for FOV and positioning / reference)
    par.dicomdir = char(dicoms(File));%F:/ICR_Studentship/data/20191205_155631_P131/scans/205_DelRec - T1 3DVane w_o fatsat thorax 1_55 mm (5_38 duration)/DICOM/*.dcm';
    %F:\ICR_Studentship\data\20191205_155631_P131\scans\205_DelRec - T1 3DVane w_o fatsat thorax 1_55 mm (5_38 duration)\DICOM

    clear a

    %%
    if ~exist([par.workspace '/par.mat'], 'file')
        fprintf('No parameter file found in working directory!\n');
        fprintf('Extracting parameters from raw / image data: \n'); tic;

        kspace = parse_list_file(par.filename);


        % loading the Dicom reference volume (slow) can be skipped by omitting the
        % last argument in the call to extract_parameters or setting it to false, 
        % but it can be a useful comparison:  
        [par, noise, refvol] = extract_parameters(par, kspace, true);
        clear kspace

        % store results for faster access in future:
        save([par.workspace '/par.mat'], 'par');
        writecfl([par.workspace '/noise'], noise);
        fid = fopen([par.workspace '/refvol.uint16'], 'w');
        fwrite(fid, refvol, 'uint16');
        [~] = fclose(fid); clear fid
        fprintf('%g seconds\n', round(toc*100)/100);
    else
        fprintf('Loading parameters: '); tic;

        load([par.workspace '/par.mat']);
        noise = readcfl([par.workspace '/noise']);

        % the reference volume is technically not needed, but nice to look at:
        fid = fopen([par.workspace '/refvol.uint16']);
        refvol = fread(fid, 'uint16=>uint16');
        [~] = fclose(fid); clear fid
        refvol = reshape(refvol, par.X_resolution, par.Y_resolution, par.Z_resolution);

        fprintf('%g seconds\n', round(toc*100)/100);
    end

    fprintf('Loading raw data: '); tic;
    % load raw data
    kdata = load_raw_data(par);
    fprintf('%g seconds\n', round(toc*100)/100);

    %%
    % perform k-space trajectory corrections:
    if ~exist([par.workspace '/knom.cfl'], 'file')
        fprintf('No trajectory information found in working directory!\n');
        fprintf('Calculating density compensation: '); tic;

        knom = (par.kx_range(1):par.kx_range(2))' / par.FOVx / 2 ... in 1/mm
            * exp(complex(0, (0:par.ky_range(2)) * par.angle_increment + par.angle_offset));

        % this calculation of the density compensation is also slow, but good. It
        % can be precalculated and the results can be loaded from disk.

        w = density_comp(knom); % the used library might not be optimal ...
        if par.useFullGPU
            w = gpuArray(w);
            knom = gpuArray(knom);
        end

        fprintf('%g seconds\n', round(toc*100)/100);
        fprintf('Performing trajectory correction: '); tic;

        % we probably should also perform a correction for gradient delays...

        [knom, w] = corr_grad_TrACR(par, knom, w, kdata);

        writecfl([par.workspace '/knom'], knom);
        writecfl([par.workspace '/w'], w);

        fprintf('%g seconds\n', round(toc*100)/100);
    else

        knom = readcfl([par.workspace '/knom']);
        w = readcfl([par.workspace '/w']);

    end


    %%
    fprintf('Estimating coil sensitivities: '); tic;

    b1_lowres = estimate_coil_sensitivities(par, knom, w, kdata, noise);

    % for upsampling:
    b1_lowres(end+1, :, :, :) = b1_lowres(end, :, :, :);
    b1_lowres(:, end+1, :, :) = b1_lowres(:, end, :, :);
    b1_lowres(:, :, end+1, :) = b1_lowres(:, :, end, :);

    kc = permute(kdata(-par.kx_range(1)+1, :, :, :), [3, 2, 4, 1]);

    % we still haven't decided the resolution for which we want to reconstruct
    % images. A natural choice would be something in the order of the acquired
    % resolution:
    nx = ceil((par.kx_range(2) + 1) / 16) * 16;
    nz = par.kz_range(2) * 2 + 1;

    % that allows us to scale the k-space accordingly:
    k = knom * par.FOVx / nx;

    % we also should probably upsample the coil sensitivity profiles: 
    [Xq, Yq, Zq] = meshgrid((1:nx)*(size(b1_lowres, 1)-1)/nx+0.5, ...
        (1:nx)*(size(b1_lowres, 2)-1)/nx+0.5, (1:nz)*(size(b1_lowres, 3)-1)/(nz-1)+0.5);

    for i = 1:size(b1_lowres, 4)
      b1(:, :, :, i) = interp3(abs(b1_lowres(:, :, :, i)),Xq,Yq,Zq,'linear').*exp(1i.*angle(interp3(b1_lowres(:, :, :, i),Xq,Yq,Zq,'linear')));
    end
    fprintf('%g seconds\n', round(toc*100)/100);  

    %save all important data to a .mat file
    par.kdata_size = size(kdata);
    par.b1_size = size(b1);
    cd(par.workspace);
    Name = strcat('data', folder, '.mat');
    save(Name, 'b1', 'kc', 'kdata', 'k', 'w', 'par');
    %save data2mat.mat b1 kc kdata k w par
    %% trying to do a FFT (no)
    %{
    load data20191205_155631_P131.mat

    kdata = fft(permute(kdata, [3 1 2 4]));
    kdata = squeeze(permute(kdata, [2 3 1 4]));
    %}
    
    %%   
    % kc is the central profiles of the stack-of-stars k-space (kx=ky=0).

    % kdata is the k-space data. It is only one slice selected from the 3D stack-of-stars
    % datasets after a FFT along the kz. In order to save recon time and
    % memory, the full stack-of-stars data is not included. Email me if you
    % want a fully 3D dataset.

    % b1 is the coil sensitivity maps of the selected slice

    % k is the radial k-space trajectory and w is the corresponding density compensation
    
    if XD

        [nz,ntviews,nc]=size(kc);
        nx=size(kdata,1);
        %nz: number of slice
        %ntviews: number of acquired spokes
        %nc: number of coil elements
        %nx: readout point of each spoke (2x oversampling included)


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Respiratory motion detection

        fprintf('Respiratory motion detection: '); tic;

        % Generate the z-projection profiles
        %ZIP: Projection profiles along the Z dimension with interpolation.
        ZIP=abs(fftshift(ifft(kc,400,1),1));

        %Normalization of each projection in each coil element
        for ii=1:nc
            for jj=1:ntviews
                maxprof=max(ZIP(:,jj,ii));
                minprof=min(ZIP(:,jj,ii));
                ZIP(:,jj,ii)=(ZIP(:,jj,ii)-minprof)./(maxprof-minprof);
            end
        end
        figure,imagesc(abs(ZIP(:,:,5))),axis image,colormap(gray), axis off

        % Perform PCA on each coil element
        close all
        kk=1;clear PCs
        for ii=1:nc
            tmp=permute(ZIP(:,:,ii),[1,3,2]);
            tmp=abs(reshape(tmp,[size(tmp,1)*size(tmp,2),ntviews])');

            covariance=cov(tmp);
            [tmp2, V] = eig(covariance);
            V = diag(V);
            [junk, rindices] = sort(-1*V);
            V = V(rindices);
            tmp2 = tmp2(:,rindices);
            PC = (tmp2' * tmp')';

            % Take the first two principal components from each coil element.
            for jj=1:2
                tmp3=smooth(PC(:,jj),6,'lowess'); % do some moving average smoothing

                %Normalize the signal for display
                tmp3=tmp3-min(tmp3(:));
                tmp3=tmp3./max(tmp3(:));
                PCs(:,kk)=tmp3;kk=kk+1;
        %         %plot the estimated signal
        %         imagesc(abs(ZIP(:,:,ii))),axis image,colormap(gray), axis off
        %         hold on
        %         plot(-tmp3(:)*100+220,'r')
        %         hold off
        %         pause
            end
        end    

        close all
        % Do coil clusting to find the respiratory motion signal
        % Function obtained from Tao Zhang (http://web.stanford.edu/~tzhang08/software.html)
        thresh = 0.95;
        [Res_Signal, cluster] = CoilClustering(PCs, thresh);

        %Normalize the signal for display
        Res_Signal=Res_Signal-min(Res_Signal(:));
        Res_Signal=Res_Signal./max(Res_Signal(:));

        % Plot the respiratory motion signal on the projection profiles
        imagesc(abs(ZIP(:,:,5))),axis image,colormap(gray), axis off,title('Respiratory Motion')
        hold on
        plot(Res_Signal(:)*100+120,'r')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
        hold off
        fprintf('%g seconds\n', round(toc*100)/100); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data sorting

        fprintf('Data sorting: '); tic;

        ntres=4;% number of respiratory phases
        nline=floor(ntviews/ntres);

        % Sort the k-space data and trajectory according to respiratory motion 
        [~,index]=sort(Res_Signal,'descend');
        kdata=kdata(:,index,:, :);
        k=k(:,index,:);
        w=w(:,index);
        kdata=kdata.*repmat(sqrt(w),[1,1,nz,nc]);

        clear kdata_u k_u w_u Res_Signal_u 
        for ii=1:ntres
            kdata_u(:,:,:,:,ii)=kdata(:,(ii-1)*nline+1:ii*nline,:, :);
            k_u(:,:,ii)=k(:,(ii-1)*nline+1:ii*nline);
            w_u(:,:,ii)=w(:,(ii-1)*nline+1:ii*nline);
        end

        fprintf('%g seconds\n', round(toc*100)/100); 
   
        %%
        % Recon

        par.noGPU = false;

        fprintf('Constructing NUFFT operator: '); tic;
        if par.useFullGPU
            %unsupported
        elseif par.noGPU
            %unsupported: use the Demo_XDGRASP_NonContrast_3D script 
        else
            param.E=MCgpuNUFFT3D(double(k_u),double(w_u),double(b1),nz);    
        end
        fprintf('%g seconds\n', round(toc*100)/100); 

        param.y = double(kdata_u);

        fprintf('Gridding 4D-MRI reconstruction: '); tic;
        recon_cs=param.E'*param.y;
        data_gridding=recon_cs/max(abs(recon_cs(:)));
        fprintf('%g seconds\n', round(toc*100)/100);

        param.TV_dim1=TV_Temp3D;
        param.TVWeight_dim1=max(abs(recon_cs(:)))*0.02;
        param.TVWeight_dim2=0;
        param.nite = 4;param.display=1;

        %clc
        tic
        for n=1:2
            recon_cs = CSL1NlCg_XDGRASP(recon_cs,param);
        end
        time=toc;
        time=time/60;
        fprintf('%g min, %g seconds\n', floor(time), round(mod(time, 60)*100)/100);
        data_xdgrasp=recon_cs/max(abs(recon_cs(:)));

        %figure, imshow3(abs(data_gridding),[0 .8],[1,4]);title('Grdding Motion State 1-4')
        %figure, imshow3(abs(data_xdgrasp),[0 .8],[1,4]);title('XD-GRASP Motion State 1-4')

        as(data_gridding)
        as(data_xdgrasp)
    
    end
    %add a log in the table using par2table
    cd(homedir)
    table = par2table(table, par);
    clearvars -except filenames dicoms XD homedir table par
end
save paramTable.mat table % this adds a table which contains all the parameters