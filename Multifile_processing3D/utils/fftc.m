function raw = fftc(raw, d)

    %raw = ifftshift(raw, d); % Philips doesn't shift in image space (!)
    raw = fft(raw, [], d);
    raw = fftshift(raw, d);

    raw = raw / sqrt(size(raw, d));

end