function cmap = estimate_coil_sensitivities(par, knom, w, kdata, varargin)

  nc = par.number_of_coil_channels;

  if nargin > 4
      rn = cov(varargin{1});
  else
      rn = eye(nc);
  end

  % the acquired resolution is given by FOV / nx, where nx is:
  nx = par.kx_range(2) + 1;

  % the number of acquired slices is given by:
  nz = par.kz_range(2) * 2 + 1;
  nky = par.ky_range(2) + 1;

  nx_low = ceil(nx / 16) * 8;
  nz_low = floor(nz / 2);

  kx_range = max(1, -par.kx_range(1) + 1 - floor(nx_low / 2)):-par.kx_range(1) + nx_low - floor(nx_low / 2);
  kz_range = max(1, -par.kz_range(1) + 1 - floor(nz_low / 2)):-par.kz_range(1) + nz_low - floor(nz_low / 2);
  nkx_low = numel(kx_range);

  k_low = knom(kx_range, :) * par.FOVx / nx_low;
  k_low = [real(k_low(:)) imag(k_low(:))]';
  w_low = w(kx_range, :);

  if par.useFullGPU
      kdata_low = complex(zeros(nkx_low, nky, nz_low, nc, 'single', 'gpuArray'));
  else
      if par.noGPU
          kdata_low = complex(zeros(nkx_low, nky, nz_low, nc, 'double'));
      else
          kdata_low = complex(zeros(nkx_low, nky, nz_low, nc, 'single'));
      end
  end
  
  kdata_low(:, :, end-numel(kz_range)+1:end, :) = kdata(kx_range, :, kz_range, :);
  kdata_low = ifftc(kdata_low, 3);
  kdata_low = reshape(kdata_low, nkx_low * nky, numel(kz_range) * nc);

  if par.useFullGPU
      NUFFT = NUFFT_GPU(reshape(complex(k_low(1, :), k_low(2, :)), size(w_low)), w_low, [nx_low nx_low]);
      x = complex(zeros(nx_low, nx_low, nz_low * nc, 'single', 'gpuArray'));
  else
      if par.noGPU
          w_low = double(w_low);
          x = complex(zeros(nx_low, nx_low, nz_low * nc, 'double'));
          NUFFT = NUFFT_irt(reshape(double(complex(k_low(1, :), k_low(2, :))), size(w_low)), double(w_low), [nx_low, nx_low]);
      else
          x = complex(zeros(nx_low, nx_low, nz_low * nc, 'single'));
          NUFFT = gpuNUFFT(k_low, w_low, 1.5, 5, 18, [nx_low, nx_low], []);
      end
  end

  res = (NUFFT * x - kdata_low .* sqrt(w_low(:)));
  g = NUFFT' * res;
  alpha = 2;
  x = x - alpha * g;
  
  x = reshape(x, size(x, 1), size(x, 2), size(x, 3) / nc, nc);
  
  %original
  [~, cmap, ~] = adapt_array_3d(x, rn);
  % alternative:
  %[~, ~, cmap] = adapt_array_3d(x, rn);
  %my fix:
  %[~, cmap] = adapt_array_3d(x, rn);
  
end
