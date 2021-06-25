function raw = load_raw_data(par)
% LOAD_RAW_DATA reads raw data from a Philips .data file starting at
% par.offset

nkx = -par.kx_range(1) + par.kx_range(2) + 1;
nky = -par.ky_range(1) + par.ky_range(2) + 1;
nkz = -par.kz_range(1) + par.kz_range(2) + 1;
nc = par.number_of_coil_channels;
raw = complex(zeros(nkx, nky, nkz, nc, 'single'));

fid = fopen(par.filename);
fseek(fid, par.offset, 'bof');

tmp = zeros(2*nkx, nc, 'single');
kz_offset = -par.kz_range(1) + 1;

for line = 1:numel(par.ky)
    
    tmp(:) = fread(fid, [2*nkx, nc], 'float32=>single');
   
    raw(:, floor((line - 1) / nkz) + 1, kz_offset + par.kz(line), :) = ...
        complex(tmp(1:2:end, :), tmp(2:2:end, :));        
    
end

if par.useFullGPU
    raw = gpuArray(raw);
end

fclose(fid);
