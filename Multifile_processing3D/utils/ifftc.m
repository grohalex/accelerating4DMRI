function raw = ifftc(raw, d)

    raw = ifftshift(raw, d);
    raw = ifft(raw, [], d);
    %raw = fftshift(raw, d); % Philips doesn't shift in image space (!)

    raw = raw * sqrt(size(raw, d));

end