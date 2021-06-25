function w = density_comp(k)

    % assumes k is complex, first dim nkx, second sim nky (angles)

    nky = size(k, 2);

    uk = unique(k(:));
    P = [real(uk) imag(uk)];
    [v,c] = voronoin(P);

    % the outmost k-space points have infinite volume associated ... 
    for j = 1:numel(c)
        x = v(c{j},:);
        if ~any(isinf(x))
          [~, volume(j)] = convhulln(x);
        else
          volume(j) = Inf;
        end
    end

    w = zeros(size(k));
    for j = 1:numel(volume)
        ind = find(k == uk(j));
        % the normalization matters for the case where no delay correction 
        % was used, because the k-space center points are identical:  
        w(ind) = volume(j) / numel(ind); 
    end

    % the k-space center point (without shift correction) is identical between
    % spokes, therefore the volume needs to be shared:
    %w_rec(-kspace_properties.kx_range(1) + 1, :) = w_rec(-kspace_properties.kx_range(1) + 1, :) / nky;

    % interpolating at the beginning and end, where the algorithm fails:
    w(1, :) = 2*w(2, :) - w(3, :);
    w(end, :) = 2*w(end-1, :) - w(end-2, :);

    % normalization on maximum might be optional
    w = w / max(w(:));

end