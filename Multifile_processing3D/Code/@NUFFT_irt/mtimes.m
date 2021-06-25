function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    %for tt=1:size(bb,4),
        %for ch=1:size(bb,3),
            ress = reshape(nufft_adj(bb .* col(a.w), a.st), a.imSize(1), a.imSize(2), size(bb, 2));
            ress = ress /sqrt(sqrt(prod(a.imSize2)) * sqrt(prod(a.imSize)));
        %end
        %ress(:,:,tt)=sum(res.*conj(a.b1),3)./sum(abs((squeeze(a.b1))).^2,3);
        %ress(:,:,tt)=sum(res.*conj(a.b1),3);
        %clear res
    %end
    %ress=ress.*size(a.w,1)*pi/2/size(a.w,2);
else
    % Cartesian image to multicoil non-Cartesian k-space
    %for tt=1:size(bb,3),
    %    for ch=1:size(a.b1,3),
    %        res=bb(:,:,tt).*a.b1(:,:,ch); 
            ress = reshape(nufft(bb,a.st),a.dataSize(1) *a.dataSize(2), size(bb, 3));
            ress = ress.*col(a.w);
            ress = ress /sqrt(sqrt(prod(a.imSize2)) * sqrt(prod(a.imSize)));
     %   end
    %end
end

