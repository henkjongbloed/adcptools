function cmap = cmaps(q)
% Returns colormap based on quantity to plot
switch q
    case "velmap"
        cmap = brewermap(20, 'RdBu');
    case "Amap"
        cmap = brewermap(20, 'Blues');
    case "phimap"
        cmap = [brewermap(20, 'RdBu') ;flipud(brewermap(20, 'RdBu'))];
    case "salmap"
        cmap = brewermap(15, 'YlOrBr');
    otherwise
        cmap = brewermap(20, 'RdBu');
end

end