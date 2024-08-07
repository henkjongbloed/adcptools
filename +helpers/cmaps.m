function cmap = cmaps(q)
% Returns colormap based on quantity to plot
n=100;
switch q
    case "velmap"
        cmap = brewermap(n, 'RdBu');
    case "Amap"
        cmap = brewermap(n, 'Blues');
    case "phimap"
        cmap = [brewermap(n, 'RdBu') ;flipud(brewermap(n, 'RdBu'))];
    case "salmap"
        cmap = brewermap(n, 'YlOrBr');
    case "fluxmap"
        cmap = brewermap(n, 'PuBuGn');
    otherwise
        cmap = brewermap(n, 'RdBu');
end

end