function ress = mtimes(a,bb)

% these we need for both adjoint and forward operation:
n_submatrix = size(a.sp, 2);
n_posperdim = a.imSize(1) / n_submatrix;

if a.adjoint
    
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    % preallocate output (this does not yet work for multiple slices at once):
    %ress = complex(zeros([size(a.b1, 1), size(a.b1, 2), size(bb, 4)], 'single', 'gpuArray'));
    % preallocate temporary variable (for coil images):
    b = complex(zeros(size(bb, 1), size(bb, 2), a.imSize(3), size(bb, 4), size(bb, 5)));  
    res = complex(zeros(size(a.b1, 1), size(a.b1, 2), size(a.b1, 3) * size(a.b1, 4), 'single', 'gpuArray'));  
    
    for tt=1:size(bb,5)
        % the loop over coil channels is slow as hell and doesn't make use
        % of the fact that we have the same NUFFT for all coils and slices:
        %for ch=1:size(bb,3),
        %    b = bb(:,:,ch,tt).*a.w(:,:,tt);
        
        b(:,:,end-size(bb,3)+1:end, :, tt) = bb(:, :, :, :, tt) .* a.w(:, :, tt);
        b(:,:, :, :, tt) = ifftc(b(:, :, :, :, tt), 3);
        
        % the operator syntax inevitable creates a copy of the data it is
        % applied to, so there is no harm in at least performing an
        % in-place multiplication with the density compensation weights:
    end
    % to perform the matrix multiplications we reshape the input data:
    b = reshape(b, size(b, 1) * size(b, 2), size(b, 3) * size(b, 4), size(b, 5));
    % after the reshaping, respiratory phases are now third dimension:
    for tt=1:size(b,3), 
        for ss = 1:n_submatrix
            for rr = 1:n_submatrix
                res(n_posperdim*(ss-1)+1:n_posperdim*ss,n_posperdim*(rr-1)+1:n_posperdim*rr,:) = ...
                    reshape(a.st{tt} * (a.sp{tt, ss, rr} .* b(:, :, tt)), n_posperdim, n_posperdim, size(b, 2));   
            end
        end
        res = res / sqrt(prod(a.imSize2));
        %res = res .* conj(a.b1);
        
        % this needs to be changed to support simultaneous application to
        % multiple slices:
        ress(:,:,tt)=sum(reshape(res, size(a.b1)) .* conj(a.b1), 4);
        %clear res
    end
    ress=ress./sum(abs((squeeze(a.b1))).^2,4);
    ress = ress * (size(a.w,1)*pi/2/size(a.w,2)); % interesting normalization ...
else
    % Cartesian image to multicoil non-Cartesian k-space
    
    % preallocate output:
    %ress = complex(zeros([size(a.w, 1) * size(a.w, 2), size(a.b1, 3), size(bb, 3)], 'single', 'gpuArray'));
    
    % preallocate temporary variable (coil images for one respiratory phase):
    %res = complex(zeros(size(a.b1), 'single', 'gpuArray'));  
    
    for tt=1:size(bb,4),
        
        % Again, we try not to loop over coils / slices:
        %for ch=1:size(a.b1,3),
        %    res=bb(:,:,tt).*a.b1(:,:,ch); 
        
        res = bb(:,:,:,tt) .* a.b1;
        res = fftc(res, 3);
        res = reshape(res, size(res, 1), size(res, 2), size(res, 3) * size(res, 4));
        
        for ss = 1:n_submatrix
            for rr = 1:n_submatrix
                
                if ss == 1 && rr == 1
                
                    ress(:,:,tt) = conj(a.sp{tt, ss, rr}) .* (a.st{tt}' * reshape( ...
                    res(n_posperdim*(ss-1)+1:n_posperdim*ss,n_posperdim*(rr-1)+1:n_posperdim*rr,:), ...
                    n_posperdim * n_posperdim, size(res, 3)));
                
                else
                    
                    ress(:,:,tt) = ress(:, :, tt) + conj(a.sp{tt, ss, rr}) .* (a.st{tt}' * reshape( ...
                    res(n_posperdim*(ss-1)+1:n_posperdim*ss,n_posperdim*(rr-1)+1:n_posperdim*rr,:), ...
                    n_posperdim * n_posperdim, size(res, 3)));
                
                end
                
            end
        end
        
        ress(:,:,tt) = ress(:, :, tt) .* col(a.w(:, :, tt));
        %end
    end
    
    ress = ress / sqrt(prod(a.imSize2));
    ress = reshape(ress, [size(a.w, 1), size(a.w, 2), a.imSize2(3), size(a.b1, 4), size(bb, 4)]);
end

