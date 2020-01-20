function plot_assoc_test
    nsample = 1000;
    n = 1000;
    ep = 0.05;
        
    Q = cell(1,3);
    
    L = abs(bsxfun(@minus, (0:n)', (0:n)));
    Px = ones(n+1, 1);
    Px = Px / sum(Px);
    
    %Exp mech
    delta = 1;
    Qygx = exp(ep/(2*delta) * (-L));
    Qygx = bsxfun(@rdivide, Qygx, sum(Qygx, 2));
    Q{1} = Qygx;
    
    %Lap mech
    Qygx = round_laplace_mech(n+1, ep);
    Q{2} = Qygx;
    
    % Ours
    Qygx = trunc_geo_mechanism(n+1, ep);
    Qzgy = opt_remap(Qygx, L, Px);
    Q{3} = Qygx * Qzgy;
    
    X = zeros(nsample, 4);
    for s=1:nsample
        r = rand()-.5;
        pxy = [.5-r .5+r; .5+r .5-r];

        px = rand();
        py = rand();
        pxy(:,1) = px * pxy(:,1);
        pxy(:,2) = (1-px) * pxy(:,2);
        pxy(1,:) = py * pxy(1,:);
        pxy(2,:) = (1-py) * pxy(2,:);
        pxy = pxy/sum(pxy(:));

        ind = randsample(4, n, true, pxy(:));
        T = [sum(ind==1), sum(ind==3);
             sum(ind==2), sum(ind==4)];
        X(s,1) = cs(T);
        for j=1:3
            Tj = zeros(2,2);
            for k=1:4
                Tj(k) = randsample(n+1,1,true,Q{j}(T(k)+1,:))-1;
            end
            X(s,j+1) = cs(Tj);
        end

    end
   
    figure;
    tls = {'Exponential','Laplace','Our Approach'};
    for j=1:3
        subplot(1,3,j)
        
        scatter(X(:,1),X(:,j+1),'b.'); hold on
        mx = max(X(:));
        plot([0 mx],[0 mx],'k-');
        
        axis equal
        axis([0 mx 0 mx])
        grid on
        xlabel('True \chi^2','FontSize',14)
        ylabel('Private \chi^2 (\epsilon = 0.05)','FontSize',14)
        set(gca,'FontSize',14);
        title(tls{j},'FontSize',14);
        fprintf('%s: %f\n', tls{j}, corr(X(:,1),X(:,j+1)).^2);
    end
    set(gcf,'color','w');
    
    function ret = cs(tbl) % chi-square statistic
        tbl = tbl+1;
        tot = sum(tbl(:));
        E = sum(tbl,2) * sum(tbl,1) / tot;
        ret = sum(sum((tbl - E).^2 ./ E));
    end
end