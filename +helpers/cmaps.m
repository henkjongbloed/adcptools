function cmap = cmaps(q)
% Returns colormap based on quantity to plot
n=100;
switch q
    case "velmap"
        cmap = helpers.brewermap(n, 'RdBu');
    case "Amap"
        cmap = helpers.brewermap(n, 'Blues');
    case "phimap"
        cmap = [helpers.brewermap(n, 'RdBu') ;flipud(helpers.brewermap(n, 'RdBu'))];
    case "salmap"
        cmap = helpers.brewermap(n, 'YlOrBr');
    case "fluxmap"
        cmap = helpers.brewermap(n, 'PiYG');
    otherwise
        cmap = helpers.brewermap(n, 'RdBu');
end

end