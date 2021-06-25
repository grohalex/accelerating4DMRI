function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    
    b = complex(zeros(size(bb, 1), size(bb, 2), a.imSize(3), size(bb, 4)));
    
    for tt=1:size(bb,5),
        
        b(:,:,end-size(bb,3)+1:end, :) = bb(:, :, :, :, tt);

        %res(:,:,ch) = reshape((a.st{tt}' * bb(:, :, tt))/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));

        %ress(:,:,:, tt)=sum(res.*conj(a.b1),3)./sum(abs((squeeze(a.b1))).^2,3);
        
        tmp(:, :, :) = a.st{tt}' * reshape(ifftc(b, 3), size(b, 1) * size(b, 2), size(b, 3) * size(b, 4));
        ress(:, :, :, tt) = sum(reshape(tmp, size(a.b1)).*conj(a.b1),4)./sum(abs((squeeze(a.b1))).^2,4);
        %ress(:, :, :, tt) = sum(reshape(tmp, size(a.b1)).*conj(a.b1),4);
    end
    
    %ress = ress /sqrt(sqrt(a.imSize2(1) * a.imSize2(2)) * sqrt(a.imSize(1) * a.imSize(2)));
    %ress = ress /sqrt(a.imSize2(1) * a.imSize2(2));
    %ress=ress.*size(a.w,1)*pi/2/size(a.w,2);
else
    % Cartesian image to multicoil non-Cartesian k-space
    for tt=1:size(bb,4),
        
        tmp(:, :, :) = reshape(bb(:, :, :, tt) .* a.b1, size(bb, 1), size(bb, 2), size(bb, 3) *size(a.b1, 4));
        
        %for ch=1:size(a.b1,3),
        %    res=bb(:,:,tt).*a.b1(:,:,ch); 
        ress(:,:,tt) = a.st{tt} * tmp;%/sqrt(prod(a.imSize2)),a.dataSize(1),a.dataSize(2));
        %end
    end
    ress = reshape(ress, size(a.w, 1), size(a.w, 2), size(bb, 3), size(a.b1, 4), size(bb, 4));

    ress = fftc(ress, 3);
    ress = ress(:, :, end-a.imSize2(3)+1:end, :, :);
    
    %ress = ress /sqrt(sqrt(a.imSize2(1) * a.imSize2(2)) * sqrt(a.imSize(1) * a.imSize(2)));
    %ress = ress /sqrt(a.imSize2(1) * a.imSize2(2));
end

