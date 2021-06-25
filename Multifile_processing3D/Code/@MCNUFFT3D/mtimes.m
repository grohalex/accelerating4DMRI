function ress = mtimes(a,bb)

% in the 3D case bb has one more dimension, i.e. the 3rd dimension is the
% slice dimension, coils are now 4th dimension and resp phases 5th.

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    
    % because of partial Fourier sampling, we need to make sure that the
    % residuals from unsampled k-space points do not contribute to the
    % gradient (and the objective function (non-adjoint case):
    bb(:, :, ~a.pf_mask, :, :) = 0;
    
    % perform fft along slice dimension
    bb = ifft(bb, [], 3);
    
    % normalization (unitary FFT/IFFT) 
    bb = sqrt(size(bb, 3)) * bb;
    
    % now the 2D NUFFT:
    for tt=1:size(bb,5)  % 4->5 
        for ch=1:size(bb,4) % 3->4
            for slc=1:size(bb, 3)
                b = bb(:,:,slc,ch,tt).*a.w(:,:,tt);
                res(:,:,slc,ch) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));
            end
        end
        ress(:,:,:,tt)=sum(res.*conj(a.b1),4)./sum(abs((squeeze(a.b1))).^2,4); % 3->4
        %clear res AW: there is no reason to clear this variable
    end
    ress=ress.*size(a.w,1)*pi/2/size(a.w,2);
else
    % Cartesian image to multicoil non-Cartesian k-space
    for tt=1:size(bb,4), % 3->4
        for ch=1:size(a.b1,4), %3->4
            for slc=1:size(bb, 3)
                res=bb(:,:,slc, tt).*a.b1(:,:,slc,ch); 
                ress(:,:,slc,ch,tt) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize2)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);       
            end
        end
    end
    
    % perform fft along slice dimension
    ress = fft(ress, [], 3);
    
    % normalization (unitary FFT/IFFT) 
    ress = ress ./ sqrt(size(bb, 3));
    
    % make sure, these points don't contribute to objective (not sampled):
    ress(:, :, ~a.pf_mask, :, :) = 0;
    
end

