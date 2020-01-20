function plot_mem_query
    n = 1000;
    
    %% Fig a: prior
    subplot(1,3,1);    
    intvl = [0 0.001 0.853;
             0.001 0.01 0.076;
             0.01 0.05 0.023;
             0.05 0.5 0.033;
             0.5 1 0.014];

    intvl(:,3) = intvl(:,3) ./ (intvl(:,2)-intvl(:,1));
    x = intvl(:,1:2)';
    x = x(:);
    y = repmat(intvl(:,3), 1, 2)';
    y = y(:);

    x = [x; 0.1];
    y = [y; 0.014+0.0037*8];

    fy = num2str(y(1));
    y(1:2)=50;

    plot(x,y,'-','LineWidth',2);
    xlabel('Allele Frequency','FontSize',14)
    ylabel('Probability Density','FontSize',14)
    axis([0 0.1 0 50])
    set(gca,'FontSize',14,'XTick',[0 0.01 0.05 0.1],'XTickLabel',{'0','0.01','0.05','>0.1'},'YTick',[0 10 20 30 40 50],'YTickLabel',{'0','10','20','30','',fy});
    grid on
    set(gcf,'color','w');
    
    %% Fig b: linear
    subplot(1,3,2);
    eta = 100;
    loss = [eta; (1:n)'];
    plot_subfigure(eta, loss)
    
    %% Fig c: error prob
    subplot(1,3,3);    
    eta = 1;
    loss = ones(n+1,1);
    plot_subfigure(eta, loss)

    function plot_subfigure(eta, loss)
        Lp = [0, loss(1); loss(2:end), zeros(n,1)];
        L = [Lp(:,1), repmat(Lp(:,2), [1, n])];
    
        ws = load('Px.mat'); % precomputed Px based on Figure (a)
        Px = ws.Px;
        Px_avg = Px;

        epvec = linspace(.2,2,50);
        Eloss = zeros(length(epvec),n+1);
        ElossExp = zeros(length(epvec),n+1);
        ElossLap = zeros(length(epvec),n+1);

        for i=1:length(epvec)
            ep=epvec(i);
            Qygx = trunc_geo_mechanism(n+1, ep);
            Qzgy = opt_remap(Qygx, L, Px);
            Qzgx = Qygx * Qzgy;
            Eloss(i,:)= sum(Qzgx .* L, 2);

            %Exp mech
            delta = eta;
            Qbgx = exp(ep/(2*delta) * (-Lp));
            Qbgx = bsxfun(@rdivide, Qbgx, sum(Qbgx, 2));
            ElossExp(i,:) = sum(Qbgx .* Lp, 2);

            %Lap mech
            Qygx = round_laplace_mech(n+1, ep);
            Qzgx = Qygx;
            ElossLap(i,:) = sum(Qzgx .* L, 2);
        end

        plot(epvec,ElossExp*Px_avg, 'k-','LineWidth', 2); hold on
        plot(epvec,ElossLap*Px_avg, 'b-','LineWidth', 2);
        plot(epvec,Eloss*Px_avg, 'r-','LineWidth', 2); hold off

        legend('Exponential','Laplace','Our Approach','FontSize',14)
        xlabel('Privacy Parameter (\epsilon)','FontSize', 14)
        ylabel('Expected Loss','FontSize',14)
        set(gca,'FontSize',14)
        grid on
    end
end
