 function [par, noise, varargout] = extract_parameters(par, kspace, varargin)
% EXTRACT_DATA uses the output from parsing the list file to find data
% dimensions and k-space trajectory.
%
% Andreas Wetscherek (a.wetscherek@icr.ac.uk 06/07/2020)

if nargin > 2
    loadimages = varargin{1};
else
    loadimages = false;
end

kspace.kspace_properties.filename = par.filename;
kspace.kspace_properties.dicomdir = par.dicomdir;
kspace.kspace_properties.workspace = par.workspace;
kspace.kspace_properties.useFullGPU = par.useFullGPU;
kspace.kspace_properties.noGPU = par.noGPU;

par = kspace.kspace_properties;

nc = par.number_of_coil_channels;
nky = par.ky_range(2) - par.ky_range(1) + 1;
nkz = par.kz_range(2) - par.kz_range(1) + 1;

ky = kspace.ky(nc+1:end); % the first nc scans are noise adjustment scans
ky = reshape(ky, nc * nkz, nky);
par.ky = ky(1:nc:end, :);
% order in which ky and kz coordinates are acquired:

if (double(ky(1,2)) / nky * 360) > 111 && (double(ky(1,2)) / nky * 360) < 112
    disp('Found golden angle trajectory');
    par.angle_increment = pi * 2 / (1 + sqrt(5)); % golden angle
    par.angle_offset = pi/2;
else
    error('unrecognized trajectory');
end

kz = kspace.kz(nc+1:end);
kz = reshape(kz, nc * nkz, nky);
par.kz = kz(1:nc:end, :);

if kz(1) == 0
    disp('Found centric order of kz partitions');
else
    error('unrecognized partition order');
end

% find offset of actual scan data:
for line = 1:numel(kspace.typ)
    if strcmp(kspace.typ{line}, 'STD') % i.e. not noise adjustment scan
        par.offset = kspace.offset(line);
        break;
    elseif ~strcmp(kspace.typ{line}, 'NOI')
        error('unrecognized scan type');
    end
end

% load noise adjustment data:
fid = fopen([par.filename(1:end-4) 'data']);
noise = fread(fid, [kspace.size(1) / 4 nc], 'float32=>single');
fclose(fid);
noise = complex(noise(1:2:end, :), noise(2:2:end, :));

disp('Reading Dicom Images');

list = dir(par.dicomdir);
for f = 1:numel(list)        
    info = dicominfo([list(f).folder '/' list(f).name]);        
    if (info.Width == par.X_resolution)            
        if ~isfield(par, 'PixelBW')
            par.PixelBW = info.PixelBandwidth;
        end
            
        if ~isfield(par, 'PixelDim')
            par.PixelDim = info.PixelSpacing;
        end
            
        if ~isfield(par, 'TE')
            par.TE = info.EchoTime;
        elseif ~ismember(par.TE, info.EchoTime)
            par.TE = [par.TE info.EchoTime];
        end
    end
    
    if ~isfield(par, 'slicepos')
        par.slicepos = info.ImagePositionPatient';
    elseif ~ismember(par.slicepos, info.ImagePositionPatient', 'rows') 
        par.slicepos = [par.slicepos; info.ImagePositionPatient'];
    end
    
    if ~loadimages
        break;
    end
    
    refvol(:, :, ismember(par.slicepos, info.ImagePositionPatient', 'rows'), ...
        par.TE == info.EchoTime) = dicomread([list(f).folder '/' list(f).name]);
end

if loadimages
    [~, I] = sort(par.slicepos(:, 3), 'descend');
    refvol = refvol(:, :, I);
    par.slicepos = par.slicepos(I, :);
end

par.FOVx = par.PixelDim(1) * par.X_resolution;
par.FOVy = par.PixelDim(2) * par.Y_resolution;

fprintf('Total time: ');

varargout{1} = refvol;
