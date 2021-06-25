
%%
% THIS PART IS JUST TO EXTRACT THE RAW DATA FROM THE PHILIPS FILES.
% IT ONLY NEEDS TO BE EXECUTED ONCE.

% this extracts all information from the raw data file pair (.list/.data)
kspace = loadRawKspacePhilips('F:\exchange\MRL_4D\20191120_125927_P131\raw_000.list');

% number of coil channels:
nc = double(max(kspace.chan)) + 1;

% number of k-space data points (ky is the radial angle basically)
nkx = kspace.kspace_properties.kx_range(2) - kspace.kspace_properties.kx_range(1) + 1;
nky = kspace.kspace_properties.ky_range(2) - kspace.kspace_properties.ky_range(1) + 1;
nkz = kspace.kspace_properties.kz_range(2) - kspace.kspace_properties.kz_range(1) + 1;

% order in which ky and kz coordinates are acquired:
ky = kspace.ky(nc+1:end); % the first nc scans are noise adjustment scans
ky = reshape(ky, nc * nkz, nky);
ky = ky(1, :);

% note that the kz loop is inside the ky loop. First all partitions (kz)
% are acquired for a particular angle, before moving to the next one.
%kz = kspace.kz(nc+1:end);
%kz = reshape(kz, nc * nkz, nky);
%kz = kz(1:nc:end, 1);

raw = complex(zeros(nkx, nc, nkz, nky, 'single'));
noise = complex(zeros(numel(kspace.complexdata{1}), nc, 'single'));

for line = 1:numel(kspace.complexdata)
    if strcmp(kspace.typ{line}, 'STD') % i.e. not noise adjustment scan
        raw(:, kspace.chan(line) + 1, ...
            kspace.kz(line) - kspace.kspace_properties.kz_range(1) + 1, ...
            kspace.ky(line) == ky) = kspace.complexdata{line};
    else
        noise(:, line) = kspace.complexdata{line};
    end
end

par = rmfield(kspace.kspace_properties, {'number_of_mixes', 'number_of_encoding_dimensions', ...
    'number_of_dynamic_scans', 'number_of_echoes', 'number_of_locations', ...
    'number_of_extra_attribute_1_values', 'number_of_extra_attribute_2_values', ...
    'number_of_signal_averages', 'number', 'X_direction', 'Y_direction', 'Z_direction', ...
    'kx_oversample_factor', 'ky_oversample_factor', 'number_of_cardiac_phases'});

writecfl('F:/alex/raw', raw);
writecfl('F:/alex/noise', noise);

par.Ne = 1; % only one echo
[refvol, par] = loadDicomVolume('F:\exchange\MRL_4D\20191120_125927_P131\scans\203_DelRec - T1 3DVane w_o fatsat thorax 1_55 mm (5_38 duration)\DICOM\*.dcm', par);

par.FTlen = round(par.Z_resolution * par.kz_oversample_factor);

par.FOVx = par.PixelDim(1) * par.X_resolution;
par.FOVy = par.PixelDim(2) * par.Y_resolution;

par.angle_increment = pi * 2 / (1 + sqrt(5)); % golden angle
par.angle_offset = pi/2;

fid = fopen('F:/alex/par.json', 'w');
fwrite(fid, jsonencode(par));
fclose(fid);

fid = fopen('F:/alex/refvol.uint16', 'w');
fwrite(fid, uint16(refvol), 'uint16');
fclose(fid);

clear

%%
% LOADING DATA

addpath('util')
addpath(genpath('../gpuNUFFT/gpuNUFFT'));
addpath(genpath('../arrShow'));
currdir = pwd;
cd('../mirt');
setup
cd(currdir);
clear 

fid = fopen('F:/alex/par.json');
par = jsondecode(fread(fid, '*char')');
[~] = fclose(fid); clear fid

fid = fopen('F:/alex/refvol.uint16');
refvol = fread(fid, '*uint16');
[~] = fclose(fid); clear fid

refvol = reshape(refvol, par.X_resolution, par.Y_resolution, par.Z_resolution);
[~, I] = sort(par.slicepos(:, 3), 'descend');
refvol = refvol(:, :, I);
par.slicepos = par.slicepos(I, :); clear I

noise = readcfl('F:/alex/noise');
raw = readcfl('F:/alex/raw');

%

[nkx, nc, nkz, nky] = size(raw);

% this is k-space position in absolute units:
knom = (par.kx_range(1):par.kx_range(2))' / par.FOVx / 2 ... in 1/mm
  * exp(complex(0, (0:nky-1) * par.angle_increment + par.angle_offset));

if exist('F:/alex/w.dbl', 'file')
    
    fid = fopen('F:/alex/w.dbl');
    w = fread(fid, 'double');
    [~] = fclose(fid); clear fid
    
else
    
    % actually acquired in-plane spatial resolution: ~1.5 mm = 500 mm / 332
Nx = par.kx_range(2) + 1;
dx = par.FOVx / Nx;

% to use it with typical NUFFT / FFT operators, we need to choose the
% resolution, that we want to reconstruct:
k = knom * dx;

% aim: correct for oversampling of k-space center
w = density_comp(k); % this can often be a simple Ram-Lak filter, 
% ... but it gets more difficult for partial echo and when trajectory
% corrections are applied.

    fid = fopen('F:/alex/w.dbl', 'w');
    fwrite(fid, w, 'double');
    [~] = fclose(fid);

end
q.data = permute(raw, [1, 4, 3, 2]); clear raw
[q.npts, q.nproj, q.npart, q.ncoils] = size(q.data);
q.dim = par.kx_range(2) + 1;
q.dimz = 2 * par.kz_range(2) + 1; %AW
q.fov = par.FOVx;    %AW field of view in mm
q.rn = cov(noise);
q.knom = knom;
q.w = w(:);
q.useGPU = 0;
clearvars -except q par

%%
% IMAGE RECONSTRUCTION USING POCS:

if q.useGPU
  % test gridding reconstruction:
  q.osf = 2; % oversampling: 2 1.5 1.25
  q.wg = 7; % kernel width: 3 5 7
  q.sw = 8; % parallel sectors' width: 8 12 16
  % ATTENTION: gpuNUFFT takes w as argument, but applies internally only sqrt(w).
  NUFFT_op = gpuNUFFT([real(q.knom(:)) imag(q.knom(:))]' * q.fov / q.dim, q.w, q.osf, q.wg, q.sw, [q.dim, q.dim], [], true);
else
  NUFFT_op = NUFFT(q.knom * q.fov / q.dim, sqrt(q.w), 1, [0, 0], [q.dim, q.dim], 2);
end
  
% actually acquired through-plane resolution: 3 mm
hybrid = complex(zeros([q.dim, q.dim, q.dimz, q.ncoils], 'single'));
offset = floor(q.dimz / 2) + par.kz_range(1);
hybrid(:, :, offset + (1:q.npart), :) = reshape(NUFFT_op' * (reshape(q.data, q.npts * q.nproj, q.npart * q.ncoils) .* sqrt(q.w)), q.dim, q.dim, q.npart, q.ncoils);

% prepare for pocs script:
hybrid = fftshift(fftshift(fft(fft(ifftshift(ifftshift(hybrid, 1), 2), [], 1), [], 2), 1), 2);
[q.imgcart, kspFull] = pocs(permute(hybrid, [4, 1, 2, 3]), 20, true);
q.imgcart = permute(q.imgcart, [2, 3, 4, 1]);

q.imgcart = ifftshift(q.imgcart, 3);

% test adaptive coil combination:

recon = complex(zeros([q.dim, q.dim, q.dimz], 'single'));
b1 = complex(zeros([q.dim, q.dim, q.dimz, q.ncoils], 'single'));

for iz = 1:q.dimz
    [recon(:, :, iz), b1(:, :, iz, :)] = adapt_array_2d(squeeze(q.imgcart(:, :, iz, :)), q.rn, true);
end; clear iz

% save another result:
iter = 0;
while (exist(['F:/alex/recon_' num2str(iter) '.cfl'], 'file'))
    iter = iter + 1;
end
writecfl(['F:/alex/recon_' num2str(iter)], recon);

%% 
% choose a slice:
slice = 59;

hybrid = complex(zeros([q.npts, q.nproj, q.dimz, q.ncoils], 'single'));
offset = floor(q.dimz / 2) + par.kz_range(1);
hybrid(:, :, offset + (1:q.npart), :) = q.data;
hybrid = ifft(ifftshift(hybrid, 3), [], 3);

% HERE WOULD ROUGHLY BE WHERE CAN EXPORT THE DATA FOR XD-GRASP
% BUT PROBABLY BETTER TO RUN THE TRAJECTORY CORRECTION FIRST 

data = reshape(hybrid(:, :, slice, :), q.npts * q.nproj, q.ncoils);
imgcart = squeeze(q.imgcart(:, :, slice, :)); % this is for TrACR code, XD-GRASP needs b1
knom = q.knom * 10;
w = q.w;
p.useGPU = q.useGPU;

clearvars -except p q par data imgcart knom w

%% update knom
cd tracr
run demo_radial_TrACR.m
cd ..

%%
% correct  k-space locations:
knom_orig = (par.kx_range(1):par.kx_range(2))' / par.FOVx / 2 ... in 1/mm
  * exp(complex(0, (0:q.nproj-1) * par.angle_increment + par.angle_offset));

if exist('ie', 'var')
    knom_orig = -knom_orig * (-1)^ie;
end

kdiff = k - [real(knom_orig(:)), imag(knom_orig(:))] * 10;
kdiff = reshape(complex(kdiff(:, 1), kdiff(:, 2)), size(knom));
kdiff = kdiff(1, :);

figure; plot(angle(knom_orig(end, :)), real(kdiff), '+')
figure; plot(angle(knom_orig(end, :)), imag(kdiff), 'r+')

kdiff_opp = interp1(angle(knom_orig(end, :)), kdiff, mod(angle(knom_orig(end, :)) + 2*pi, 2*pi) - pi);
todo = isnan(kdiff_opp);
kdiff_opp(todo) = interp1(mod(angle(knom_orig(end, :)) + 2 * pi, 2*pi), kdiff, mod(angle(knom_orig(end, todo)) + pi, 2*pi));

%figure; plot(angle(knom_orig(end, :)), real(kdiff_opp), 'c+')
%figure; plot(angle(knom_orig(end, :)), imag(kdiff_opp), 'k+')

kcorr = (kdiff - kdiff_opp) / 2;

figure; plot(angle(knom_orig(end, :)), real(kcorr), 'c+')
figure; plot(angle(knom_orig(end, :)), imag(kcorr), 'k+')

kcorr = repmat(kcorr, size(knom,1), 1);
k_final = knom_orig * 10 + kcorr;

k_final = [real(k_final(:)) imag(k_final(:))];

% update nominal k-space locations
q.knom = reshape(complex(k_final(:, 1), k_final(:, 2)) / 10, size(q.knom));

% and density compensation:
q.w = density_comp(double(q.knom));
q.w = q.w(:);
clearvars -except q par

% then repeat IMAGE RECONSTRUCTION USING POCS (line 145)

%%
% compare two recon iterations side by side:
left = 0;
right = 1;

comp = readcfl(['F:/alex/recon_' num2str(left)]);
comp = repmat(comp, [1, 2, 1]);
comp(:, end/2+1:end, :) = readcfl(['F:/alex/recon_' num2str(right)]);

as(comp)

%%
% IGNORE REST OF THIS SCRIPT





% do it in 3D:

data = reshape(q.data, q.npts * q.nproj, q.npart, q.ncoils);
imgcart = q.imgcart;
knom = q.knom * 10;
w = q.w;

run demo_radial_TrACR_3D.m

q.knom = reshape(complex(k_final(:, 1), k_final(:, 2)) / 10, size(q.knom));
q.w = density_comp(double(q.knom));
q.w = q.w(:);
clearvars -except q par


