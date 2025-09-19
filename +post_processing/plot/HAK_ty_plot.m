function treal = HAK_ty_plot(D, AF)
cm = ["velmap", "salmap", "fluxmap"];

di = 1;
bi = 1;
ri=1;
treal{di} = datetime(D{di}{bi}.X.t, 'ConvertFrom', 'datenum');
tit = ["Flow velocity $[u](t,y)$ [m/s]", "Salinity $[s](t,y)$ [psu]",  "Flux $[f](t,y)$ [psu m/s]"];

[fig, TL] = prep_fig([9, 12], [3,1]);
TL.Padding = "compact";
TL.TileSpacing = "compact";
c=1;
for vi = 1:3 %flow, salt, flux
    ax = nexttile;
    dat = AF{di}{bi}{vi}{ri, c}{7}'; %6 : depth averaged

    contourf(squeeze(D{di}{bi}.X.T(:,:,1))', squeeze(D{di}{bi}.X.Y(:,:,1))', squeeze(dat), 50, 'EdgeColor','none');
    shading interp;
    hold on
    contour(squeeze(D{di}{bi}.X.T(:,:,1))', squeeze(D{di}{bi}.X.Y(:,:,1))', squeeze(dat), 5, 'LineWidth',.3, 'EdgeColor','k' );
%     plot(D{di}{bi}.X.t, D{di}{bi}.wl,'Color','k', 'LineWidth',2)
    %         shading interp;
    %set(c, 'edgecolor','none');
    if vi==1
        clim(ax, [-max(abs(dat), [], "all"), max(abs(dat), [], "all")])
        colormap(ax, flipud(helpers.cmaps('velmap')));
        ax.XTickLabel = [];
    elseif vi==2
        clim(ax, [0, max(abs(dat), [], "all")])
        colormap(ax, helpers.cmaps('salmap'));
        
    else
        clim(ax, [-max(abs(dat), [], "all"), max(abs(dat), [], "all")])
        colormap(ax, flipud(helpers.cmaps('fluxmap')));
    end
    colorbar
    title(tit(vi), 'Interpreter','latex')
        ylabel("y [m]")
    ax.FontSize = 12;
%     ax.YTickLabel = [];
    ax.XGrid = 'on';
    ax.XTickLabel = [];
end
fontname(TL, "Book Antiqua")
sgtitle("Hartel Canal: Depth-averaged", "FontName", "Book Antiqua")
% 
% 
% 
%
end