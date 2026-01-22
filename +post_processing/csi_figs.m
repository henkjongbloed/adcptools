function csi_figs(di, bidx, D, AF)
% Selection of data
tit = {"Junction S", "Junction N"};
ls = {'-', '--', '-.'};
for bi = bidx
    A{bi} = D{di}{bi}.A/1000;
end
sgtitle(tit{di}, "FontName","Book Antiqua")
c = get_c();
colororder({'k','k'})
% TL = tiledlayout(4, 1);
treal{di} = datetime(D{di}{bidx(1)}.X.t, 'ConvertFrom', 'datenum');
opts.ylim = 0;
ax{1} = prep_ax();

axis tight
ax{1}.XTickLabel = [];
ylabel("m", "FontName","Book Antiqua", "Interpreter","latex")
yyaxis left
hold on
for bi = 1 % Only one waterlevel
    plot(treal{di}, D{di}{bi}.wl,'Color','k', 'LineWidth', 2)
end
yyaxis(ax{1}, 'right')
opts.alternating_y = 0;
for bi = bidx
    plot(treal{di}, A{bi}, 'Color', c{bi}, 'LineStyle', ls{2}, 'LineWidth', 2)
end
hold off
ylabel("10$^3$ m$^2$", "FontName","Book Antiqua", "Interpreter","latex")

%     plot(treal, eta{1,di}, c{1}, 'LineWidth',1)
title("Waterlevel $\eta$/ Wetted area A", "FontName","Book Antiqua", "Interpreter", "latex")
%legend(["Waterlevel", tak_name{di}{:}],'AutoUpdate','off')
opts.dashedleft = 0;
opts.dashedright = 1;
ax{1}.XTickLabel = [];
ax{1} = fix_ax(ax{1}, opts);

%flow
ax{2}= prep_ax();hold on
for bi = bidx
    plot(treal{di} , AF{di}{bi}{1}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{1})
end
ylabel("m/s", "FontName","Book Antiqua", "Interpreter","latex")

yyaxis right
for bi = bidx
    plot(treal{di} , A{bi}.*AF{di}{bi}{1}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{2}, 'LineWidth', 2)
end
hold off
ax{2}.XTickLabel = [];
ylabel("10$^3$ m$^3$/s", "FontName","Book Antiqua", "Interpreter","latex")
%legend(tak_name{di}{:})
axis tight
title("Flow velocity $u$ (avg / int)", "FontName","Book Antiqua", "Interpreter", "latex")
ax{2} = fix_ax(ax{2}, opts);
%sumu = D{1,di}.A.* AF{1,di}{1,2} + D{2,di}.A.* AF{2,di}{1,2}- D{3,di}.A.* AF{3,di}{1,2};
%plot(taxis , sumu, 'Color', c{2}, 'LineStyle', ls{bi})

%salt
% fix_ax(ax{2}, opts)
ax{3}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, AF{di}{bi}{2}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{1}, 'LineWidth', 2)
end
ylabel("psu", "FontName","Book Antiqua", "Interpreter","latex")

yyaxis right
for bi = bidx
    plot(treal{di}, A{bi}.*AF{di}{bi}{2}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{2}, 'LineWidth', 2)
end
hold off
ylabel("10$^3$ psu m$^2$", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("Salinity $s$ (avg / int)", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = "pos";
ax{3}.XTickLabel = [];

ax{3} = fix_ax(ax{3}, opts);

%flux
% fix_ax(ax{2}, opts)
ax{4}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, AF{di}{bi}{3}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{1}, 'LineWidth', 2)
end
ylabel("psu m/s", "FontName","Book Antiqua", "Interpreter","latex")

yyaxis right
for bi = bidx
    plot(treal{di}, A{bi}.*AF{di}{bi}{3}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{2}, 'LineWidth', 2)
end
hold off
ylabel("10$^3$ psu m$^3$/s", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("Flux $f = us$ (avg / int)", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim =0;
% opts.XTick = 0;
ax{4} = fix_ax(ax{4}, opts);

end
