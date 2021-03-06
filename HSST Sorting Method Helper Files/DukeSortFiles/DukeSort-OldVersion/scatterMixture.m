function h = scatterMixture(data, z)
    colors = ['m','b','r','k','c','g','y'];
    symbols = ['.','o','x','+','*','s','d','v','^','<','>','p','h'];
    h = [];
    
    for i = min(z):max(z)
        hold on;
        format = strcat(colors(1+rem(i,numel(colors))),symbols(1+rem(3*i,numel(symbols))));
        h(i) = scatter(data(z==i,1),data(z==i,2),format);
    end
    legend;
end