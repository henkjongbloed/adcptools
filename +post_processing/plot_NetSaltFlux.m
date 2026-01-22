function plot_NetSaltFlux(didx, bidx, SF, bname)

tl = tiledlayout(2,3, "TileSpacing","compact", Padding="compact"); % One fig per dataset % two cols: SF and UU,
xt = ["$u_0s_0$", "$\overline{u_1s_1}^t$", "$\overline{u_2s_2}^y$", "$\overline{u_3s_3}^\sigma$",...
    "$\overline{u_4s_4}$", "$\overline{u_5s_5}$", "$\overline{u_6s_6}$", "$\overline{u_7s_7}$"];
yabs=0;
c = get_c();
for di = didx
    for bi = bidx
        dat{di, bi} = [SF{di}{bi}{1,1}{1,:,3}]';
        yabs = max(yabs, max(abs(dat{di, bi})));
    end
    ylims{di} = [-yabs, yabs];
end

for di = didx
    for bi = bidx
        ax=nexttile();
        fontname(ax, "Book Antiqua")

        

        xticks(0:7)
        hold on
        grid on
            b=bar(0:7, dat{di, bi}, "FaceColor", c{bi});
%         b.CData = c(bi);
        axis tight
        ylim(ylims{di})
        title([bname{di}{bi}], 'Interpreter','latex')

        if bi==1
            ylabel("psu m/s", 'Interpreter','latex')
        else
            ax.YTickLabel = [];
        end
        ax.FontSize = 10;
        if di == 2
            ax.XTickLabel = xt;
            ax.XAxis.FontSize = 12;
            
%             ax.Tick
        else
            ax.XTickLabel = [];
        end
        ax.TickLabelInterpreter = 'latex';
        
        
    end
end
fontname(tl, "Book Antiqua")
end