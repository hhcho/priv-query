function Pzgy = opt_remap(Pygx, L, Px)
    if ~exist('Px','var')
        Px = ones(size(Pygx, 1), 1) / size(Pygx, 1);
    end
    
    [n, m] = size(Pygx);
    
    Pyx = bsxfun(@times, Px, Pygx);    
    Pxgy = bsxfun(@rdivide, Pyx, sum(Pyx, 1));
    [~, opt_inds] = min(Pxgy' * L,[],2);
    Pzgy = sparse((1:length(opt_inds))', opt_inds(:), 1, n, m);
   
end