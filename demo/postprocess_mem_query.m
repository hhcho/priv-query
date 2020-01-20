% [Input]
% y: perturbed query result received from database (between 0 and n),
%    which is ep-differentially private
% n: database size
% ep: privacy parameter epsilon
% px: prior over true count x (length n + 1 vector),
%     OK to be unnormalized
% loss: loss function given as a vector of length n + 1,
%       where the i-th element corresponds to the loss incurred
%       when flipped answer is returned when true count is i
% [Output]
% z: loss minimizing guess for the query result
function z = postprocess_mem_query(y, n, ep, px, loss)
    if length(px) ~= n + 1
        error('Length of px does not equal n+1');
    end
    if length(loss) ~= n + 1
        error('Length of loss does not equal n+1');
    end
    if y < 0 || y > n
        error('Invalid value of y');
    end
    
    alpha = exp(-ep);
    post_px = px(:) .* (alpha .^ abs((0:n)' - y)); % posterior over x
    
    z = double(post_px(1) * loss(1) < dot(post_px(2:end), loss(2:end)));
end