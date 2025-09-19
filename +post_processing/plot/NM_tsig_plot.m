function NM_tsig_plot(treal, D, AF)
cm = ["velmap", "salmap", "fluxmap"];

di = 2;
bi = 1;
ri=1;
c=1;
treal{di} = datetime(D{di}{bi}.X.t, 'ConvertFrom', 'datenum');
tit = ["Flow velocity $\underline{u}(t,\sigma)$ [m/s]", "Salinity $\underline{s}(t,\sigma)$ [psu]",  "Flux $\underline{f}(t,\sigma)$ [psu m/s]"];
%            xticklabels(["$u_0$", "$\bar{u}_t^t$", "$\bar{u}_y^y$", "$\bar{u}_\sigma^\sigma$",...
%                "$\widehat{u}_{y\sigma}^{y\sigma}$", "$\underline{u}_{t\sigma}^{t\sigma}$", "$[u]_{ty}^{ty}$", "$u_{t y\sigma}$"]);
%            yticklabels(["$s_0$", "$\bar{s}_t^t$", "$\bar{s}_y^y$", "$\bar{s}_\sigma^\sigma$",...
%               "$\widehat{s}_{y\sigma}^{y\sigma}$", "$\underline{s}_{t\sigma}^{t\sigma}$", "$[s]_{ty}^{ty}$", "$s_{t y\sigma}$"]);
% leg = {["HC", "OMS", "OMN"], ["NM", "OM", "NWW"]};
residx_ = 0:7; % change this to 0:7 for the full flow field
residx = residx_ + 1;
sgtit = ["Junction S: Inst. component norms", "Junction N: Inst. component norms"];
%             for vi = 1:3 %flow, salt, flux
Zplot = D{di}{bi}.wi(D{di}{bi}.X.Z)/D{di}{bi}.B;
% Zplot2 = D{di}{bi}.la(D{di}{bi}.X.Z);
[fig, TL] = prep_fig([9, 12], [3,1]);
TL.Padding = "compact";
TL.TileSpacing = "compact";

for vi = 1:3 %flow, salt, flux
    ax = nexttile;
    dat = AF{di}{bi}{vi}{ri, c}{6}; %6 : width averaged

    contourf(squeeze(D{di}{bi}.X.T(:,1,:)), squeeze(Zplot), squeeze(dat), 50, 'EdgeColor','none');
    shading interp;
    hold on
    contour(squeeze(D{di}{bi}.X.T(:,1,:)), squeeze(Zplot), squeeze(dat), 5, 'LineWidth',.3, 'EdgeColor','k' );
    plot(D{di}{bi}.X.t, D{di}{bi}.wl,'Color','k', 'LineWidth',2)
    %         shading interp;
    %set(c, 'edgecolor','none');
    if vi==1
        clim(ax, [-max(abs(dat), [], "all"), max(abs(dat), [], "all")])
        colormap(ax, flipud(helpers.cmaps('velmap')));
    elseif vi==2
        clim(ax, [0, max(abs(dat), [], "all")])
        colormap(ax, helpers.cmaps('salmap'));
    else
        clim(ax, [-max(abs(dat), [], "all"), max(abs(dat), [], "all")])
        colormap(ax, flipud(helpers.cmaps('fluxmap')));
    end
    colorbar
    title(tit(vi), 'Interpreter','latex')
        ylabel("z [m]")
    ax.FontSize = 12;
    ax.XTickLabel = [];
    ax.XGrid = 'on';
end
fontname(TL, "Book Antiqua")
sgtitle("New Meuse: Width-averaged", "FontName", "Book Antiqua")
end