function plot_NetSaltBalance(didx, bidx, SF, bname)

tl = tiledlayout(2,1, "TileSpacing","compact", Padding="compact"); % One fig per dataset % two cols: SF and UU,
xt = ["$u_0s_0$", "$\overline{u_1s_1}^t$", "$\overline{u_2s_2}^y$", "$\overline{u_3s_3}^\sigma$",...
    "$\overline{u_4s_4}$", "$\overline{u_5s_5}$", "$\overline{u_6s_6}$", "$\overline{u_7s_7}$", "$f_0$"];
yabs=0;
c = get_c();
bsign = {[1, -1, 1], [-1, -1, 1]};
for di = didx
    for bi = bidx
        dat{di}(bi,1:8) = bsign{di}(bi).*[SF{di}{bi}{1,1}{1,:,3}];
        dat{di}(bi,9) = sum(dat{di}(bi,1:8)); % Sum of mechanisms per branch
        %         yabs = max(yabs, max(abs(dat{di, bi})));
    end
    %     ylims{di} = [-yabs, yabs];
end

for di = didx
    summ{di}  = sum(dat{di}, 1);
    %sumb{di} =
    % sumt{di} = sum(summ{di});
end
tit = {"Junction S", "Junction N"};
for di = flip(didx)
    ax=nexttile();
    title(tit{di})
    fontname(ax, "Book Antiqua")



    xticks(0:8)
    hold on
    grid on
    b=bar(0:8, dat{di}, "stacked");
    bs = bar(0:8, summ{di}, .5, 'FaceColor', 'k', 'FaceAlpha', 1);
    axis tight

    ylabel("psu m/s", 'Interpreter','latex')
    %ax.FontSize = 10;
    ax.XTickLabel = xt;
    %ax.XAxis.FontSize = 12;
    
    ax.TickLabelInterpreter = 'latex';
    legend({bname{di}{:}, 'Sum'})
    sgtitle("Net salt balance")
end
fontname(tl, "Book Antiqua")
end