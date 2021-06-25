function [knom, w] = corr_grad_TrACR(par, knom, w, kdata)

  % the acquired resolution is given by FOV / nx, where nx is:
  nx = par.kx_range(2) + 1;

  % the number of acquired slices is given by:
  nz = par.kz_range(2) * 2 + 1;
  nc = par.number_of_coil_channels;
  nky = par.ky_range(2) + 1;

  % try to perform trajectory correction based on a low resolution image:
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
      [px, py] = ndgrid(gpuArray(-0.5:1/nx_low:0.5-0.5/nx_low));
  else
      [px, py] = ndgrid(-0.5:1/nx_low:0.5-0.5/nx_low);
      if par.noGPU
          x = complex(zeros(nx_low, nx_low, nz_low * nc, 'double'));
          NUFFT = NUFFT_irt(reshape(double(complex(k_low(1, :), k_low(2, :))), size(w_low)), double(w_low), [nx_low, nx_low]);
      else
          x = complex(zeros(nx_low, nx_low, nz_low * nc, 'single'));
          NUFFT = gpuNUFFT(k_low, w_low, 1.5, 5, 18, [nx_low, nx_low], []);
      end
  end

  res = (NUFFT * x - kdata_low .* sqrt(w_low(:)));
  fprintf('iter 0: %g\n', real(res(:)' * col(res ./ w_low(:))) / 2);

  dk = k_low; % ensures we have the same data type ...
  dk(:) = 0;

  for f = 1:2

    g = NUFFT' * res;
    alpha = 2;
    x = x - alpha * g;
    res = (NUFFT * x - kdata_low .* sqrt(w_low(:)));
    fprintf('image update %d: alpha = %g, f = %g\n', f, alpha, real(res(:)' * col(res ./ w_low(:))) / 2);

    %
    g = [sum(real(conj(1i*2*pi*(NUFFT * (x.*px))).*res), 2).'; ...
         sum(real(conj(1i*2*pi*(NUFFT * (x.*py))).*res), 2).'];
 
    g = repmat(mean(reshape(g, 2, nkx_low, nky), 2), [1, nkx_low, 1])/nx_low;

    alpha = 1/max(sqrt(sum(g(:, :).^2, 1)));

    g = reshape(g, size(k_low));
 
    %
    beta = 1e-4;
    q = 0.5;
    maxit = 30;

    y0 = res(:)' * res(:)/2;
    g0 = (g(:)' * g(:));
    
    if par.useFullGPU
      k_tmp = k_low + dk + alpha * g;
      NUFFT = NUFFT_GPU(reshape(complex(k_tmp(1, :), k_tmp(2, :)), size(w_low)), w_low, [nx_low nx_low]);
    elseif par.noGPU
      k_tmp = k_low + dk + alpha * g;
      NUFFT = NUFFT_irt(reshape(complex(k_tmp(1, :), k_tmp(2, :)), size(w_low)), w_low, [nx_low nx_low]);
    else
      NUFFT = gpuNUFFT(k_low + dk + alpha * g, w_low, 1.5, 5, 18, [nx_low, nx_low], []);
    end

    res = (NUFFT * x - kdata_low .* sqrt(w_low(:)));
    y1 = res(:)' * res(:)/2;
    it = 0;

    while y1 > y0 - beta * alpha * g0 && it < maxit
        it = it + 1;
        alpha = alpha * q;
        if par.useFullGPU
            k_tmp = k_low + dk + alpha * g;
            NUFFT = NUFFT_GPU(reshape(complex(k_tmp(1, :), k_tmp(2, :)), size(w_low)), w_low, [nx_low nx_low]);
        elseif par.noGPU
            k_tmp = k_low + dk + alpha * g;
            NUFFT = NUFFT_irt(reshape(complex(k_tmp(1, :), k_tmp(2, :)), size(w_low)), w_low, [nx_low nx_low]);
        else
            NUFFT = gpuNUFFT(k_low + dk + alpha * g, w_low, 1.5, 5, 18, [nx_low, nx_low], []);
        end
        res = (NUFFT * x - kdata_low .* sqrt(w_low(:)));
        y1 = res(:)' * res(:)/2;
    end

    fprintf('trajectory update %d: alpha = %g, f = %g\n', f, alpha, real(res(:)' * col(res ./ w_low(:))) / 2);

    dk = dk + alpha * g;

  end

  traj_corr = reshape(complex(dk(1, :), dk(2, :)), nkx_low, nky);
  traj_corr = traj_corr(1, :);

  angles = atan2(imag(knom(1, :)), real(knom(1, :)));

  [asort, I] = sort(angles, 'ascend');
  fw = 39;
  ex = fw/2-0.5;
  traj_corr_ma(I) = movmean([traj_corr(I(end-ex+1:end)) traj_corr(I) traj_corr(I(1:ex))], fw, 'Endpoints','discard');

  figure; plot(angles, real(traj_corr), '+'); 
  hold on; plot(angles, real(traj_corr_ma), '+'); 

  figure; plot(angles, imag(traj_corr), '+'); 
  hold on; plot(angles, imag(traj_corr_ma), '+'); 
    
  [a_opp, J] = sort(mod(angles + 2 * pi, 2 * pi) - pi, 'ascend');
  traj_opp = traj_corr_ma(J);

  traj_int = interp1([a_opp(end)-2*pi a_opp a_opp(1)+2*pi], [traj_opp(end) traj_opp traj_opp(1)], asort);
  traj_opp(I) = traj_int;

  plot(angles, imag(traj_corr_ma), '+');
  hold on
  plot(angles, imag(-traj_opp), '+');

  traj_corr_final = (traj_corr_ma - traj_opp) / 2;

  traj_corr_final = repmat(traj_corr_final, size(knom, 1), 1);

  knom = knom + double(traj_corr_final * nx_low / par.FOVx);
  
  if par.useFullGPU
      w = gpuArray(density_comp(gather(knom)));
  else
      w = density_comp(knom);
  end

end
