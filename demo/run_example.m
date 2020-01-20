function run_example
    %% Cohort discovery
    fprintf('## Cohort discovery ##\n')
    fprintf('Database size: 1000000, True count result: 1000, Epsilon: 0.1\n');
    n = 1e6;
    x = 1000;
    ep = 0.1;
    
    tic
    y = trunc_geo_mechanism(x, n, ep);
    fprintf('\t[Database->User] Perturbed count: %d, ', y);
    toc
    
    fprintf('Prior distribution: Uniform, Loss function: Asymmetric linear\n');
    px = ones(1, n+1);
    loss_fn = @(delta) ((delta>=0) .* (2 * delta) + (delta<0) .* (-delta));
    
    tic
    z = postprocess_cnt_query(y, n, ep, px, loss_fn);
    fprintf('\t[User postprocessing] Optimal query result: %d, ', z);
    toc
    
    %% Variant lookup
    fprintf('\n## Variant lookup ##\n')
    fprintf('Database size: 1000000, True count result: 50, Epsilon: 0.1\n');
    n = 1e6;
    x = 50;
    ep = 0.1;
    
    tic
    y = trunc_geo_mechanism(x, n, ep);
    fprintf('\t[Database->User] Perturbed count: %d, ', y);
    toc
    
    fprintf('Prior distribution: Uniform, Loss function: Linear\n');
    px = ones(1, n+1);
    loss = [1 1:n];
    
    tic
    z = postprocess_mem_query(y, n, ep, px, loss);
    ans = {'no', 'yes'};
    fprintf('\t[User postprocessing] Optimal query result: %d (%s), ', z, ans{z+1});
    toc

end