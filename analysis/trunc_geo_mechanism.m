function Pygx = trunc_geo_mechanism(n, ep)

    alpha = exp(-ep);
    norm_factor = (1-alpha) / (1+alpha);
    Pygx = norm_factor * (alpha .^ abs(bsxfun(@minus, (1:n)', 1:n)));
    
    % Truncated parts
    accum = (1 / (1 + alpha)) * (alpha .^ (1:n)');
    Pygx(:,1) = Pygx(:,1) + accum;
    Pygx(:,end) = Pygx(:,end) + accum(end:-1:1);
        
end