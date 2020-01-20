function Qygx = round_laplace_mech(n, ep)
    r_yminusx = bsxfun(@plus, -(1:n)', (1:n) + .5);
    l_yminusx = bsxfun(@plus, -(1:n)', (1:n) - .5);
    r_yminusx(:,end) = inf;
    l_yminusx(:,1) = -inf;
    
    sgn = r_yminusx > 0;
    rF = sgn - (sgn*2-1) .* 0.5 .* exp(-ep*abs(r_yminusx));
    
    sgn = l_yminusx > 0;
    lF = sgn - (sgn*2-1) .* 0.5 .* exp(-ep*abs(l_yminusx));
    
    Qygx = rF - lF;
    Qygx = bsxfun(@rdivide, Qygx, sum(Qygx, 2));
end