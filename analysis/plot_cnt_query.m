function plot_cnt_query

    n = 1000;    
    alpha = [0.5 0.5];
    beta = [3 1];
    
    cost_fn = @(d) beta(1) * (d >= 0) .* (d .^ alpha(1)) + beta(2) * (d < 0) .* ((-d) .^ alpha(2));
    L = cost_fn(bsxfun(@minus, 1:(n+1), (1:(n+1))'));
    
    %% Figure a: loss
    figure
    subplot(1,3,1)
    x = linspace(0,100,2000)-50;
    xc = cost_fn(x);
    
    plot(x,xc,'-','Linewidth',2)
    xlabel('Absolute error','FontSize',14)
    ylabel('Loss','FontSize',14)
    set(gca,'FontSize',14)
    grid on
    
    %% Figure b: exp loss
    subplot(1,3,2)
    Px = ones(n+1, 1);
    Px = Px / sum(Px);    
    Px_avg = Px;
    
    epvec = linspace(.2,2,50);
    
    Eloss = zeros(length(epvec),n+1);
    ElossExp = zeros(length(epvec),n+1);
    ElossLap = zeros(length(epvec),n+1);
    
    Qzgx_cache = cell(1,3);
    
    for i=1:length(epvec)
        ep = epvec(i);
        
        Qygx = trunc_geo_mechanism(n+1, ep);
        Qzgy = opt_remap(Qygx, L, Px);
        Qzgx = Qygx * Qzgy;
        Eloss(i,:)= sum(Qzgx .* L, 2);
        if i==1; Qzgx_cache{3} = Qzgx; end

        %Exp mech
        delta = max(beta);
        Qygx = exp(ep/(2*delta) * (-L));
        Qygx = bsxfun(@rdivide, Qygx, sum(Qygx, 2));
        ElossExp(i,:) = sum(Qygx .* L, 2);
        if i==1; Qzgx_cache{1} = Qygx; end

        %Lap mech
        Qygx = round_laplace_mech(n+1, ep);
        Qzgx = Qygx;
        ElossLap(i,:) = sum(Qzgx .* L, 2);
        if i==1; Qzgx_cache{2} = Qzgx; end
    end
    
    plot(epvec,min(5,ElossExp*Px_avg), 'k-','LineWidth', 2); hold on
    plot(epvec,ElossLap*Px_avg, 'b-','LineWidth', 2);
    plot(epvec,Eloss*Px_avg, 'r-','LineWidth', 2); hold off
    legend('Exponential','Laplace','Our Approach','FontSize',14)
    xlabel('Privacy Parameter (\epsilon)','FontSize', 14)
    ylabel('Expected Loss','FontSize',14)
    set(gca,'FontSize',14,'YTick',0:5,'YTickLabel',{'0','1','2','3','4','>5'})
    axis([.2 2, 0 5])
    grid on
    
    %% Figure c: output prob
    subplot(1,3,3)
    
    cl = {'k-','b-','r-'};
    for i=1:3
       plot(0:100,Qzgx_cache{i}(51,1:101), cl{i},'LineWidth', 2); hold on
    end
    legend('Exponential','Laplace','Our Approach','FontSize',14)
    xlabel('Released Count','FontSize', 14)
    ylabel('Probability','FontSize',14)
    set(gca,'FontSize',14);
    grid on
    set(gcf,'color','w');

end