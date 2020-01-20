% [Input]
% x: true count
% n: database size
% ep: privacy parameter epsilon (non-negative)
% [Output]
% y: ep-differentially private release of x
%    using truncated geometric mechanism
function y = trunc_geo_mechanism(x, n, ep)
    if x < 0 || x > n
        error('Invalid value of x');
    end
    if ep < 0
        error('Invalid value of ep');
    end
    
    if ep < eps
        y = randi(n + 1) - 1;
        return
    end
    
    alpha = exp(-ep);
    
    % flip a coin to determine membership in
    % [0, inf) or (-inf, -1], both of which
    % are sampled using geometric distribution
    ispos = rand() < 1 / (1 + alpha); 
    
    delta = geornd(1 - alpha);
    if ~ispos
        delta = - delta - 1;
    end
    y = max(min(x + delta, n), 0);

end