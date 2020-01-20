% [Input]
% y: perturbed query result received from database (between 0 and n),
%    which is ep-differentially private
% n: database size
% ep: privacy parameter epsilon
% px: prior over true count x (length n + 1 vector),
%     OK to be unnormalized
% loss_fn: loss function given as a function handle called as
%          loss_fn(delta), where delta = released count - true count;
%          should also perform elementwise calc for a vector-valued delta
% [Output]
% z: loss minimizing guess for the query result
function z = postprocess_cnt_query(y, n, ep, px, loss_fn)
    if length(px) ~= n + 1
        error('Length of px does not equal n+1');
    end
    if y < 0 || y > n
        error('Invalid value of y');
    end
    
    alpha = exp(-ep);
    post_px = px(:) .* (alpha .^ abs((0:n)' - y)); % posterior over x
    
    L = loss_fn(n:-1:-n);
    exp_loss = fconv(L(:), post_px(:));
    [~, z] = min(exp_loss(n+1:end-n));
    z = z - 1;
end