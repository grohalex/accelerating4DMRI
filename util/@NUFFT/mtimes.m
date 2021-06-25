function ress = mtimes(a,bb)
% performs the normal nufft

if numel(size(bb)) == 2
  bb = reshape(bb, size(bb, 1), 1, size(bb, 2));
end

for n=1:size(bb,3)
b = bb(:,:,n);

if a.adjoint
	b = b(:).*a.w(:);
	res = nufft_adj(b, a.st)/sqrt(prod(a.imSize));
	res = reshape(res, a.imSize(1), a.imSize(2));
	res = res.*conj(a.phase);
	if a.mode==1
		res = real(res);
    end
    ress(:,:,n) = res;
else
	b = reshape(b,a.imSize(1),a.imSize(2));
	if a.mode==1
		b = real(b);
	end
	b = b.*a.phase;
	res = nufft(b, a.st)/sqrt(prod(a.imSize));
	res = reshape(res,a.dataSize(1),a.dataSize(2));
    ress(:,n) = res;
end

end

