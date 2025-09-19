function strat_figs(di, bidx, D, laIFB, laIFS, daIFP, daIFN)


ls = {'-', '--', '-.'};
tit = {"Junction S", "Junction N"};
c = get_c();
colororder({'k','k'})
TL = tiledlayout(3,2);
treal{di} = datetime(D{di}{bidx(1)}.X.t, 'ConvertFrom', 'datenum');
opts.ylim = 0;
fa = [.8, .7, .6];

%velocity stratification
% fix_ax(ax{2}, opts)
ax{1}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, laIFS{di}{bi}{1}{1}, 'Color', 'k', 'LineStyle', ls{1}, 'LineWidth',1.5)
    plot(treal{di}, laIFB{di}{bi}{1}{1}, 'Color', 'k', 'LineStyle', ls{2}, 'LineWidth',1.5)
    fill_between(treal{di}', laIFS{di}{bi}{1}{1}, laIFB{di}{bi}{1}{1} ,  c{bi}, fa(bi));
end
ylabel("m/s", "FontName","Book Antiqua", "Interpreter","latex")

% yyaxis right
% for bi = bidx
%     plot(treal{di}, A{bi}.*AF{di}{bi}{3}{1}{2}, 'Color', c{bi}, 'LineStyle', ls{2}, 'LineWidth',1)
% end
% hold off
% ylabel("10$^3$ psu m$^3$/s", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("$u$: $\sigma$-shear", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = 'sym';
% % opts.XTick = 0;
ax{1} = fix_ax(ax{1}, opts);
% fix_ax(ax{2}, opts)
ax{1}.XTick = [];

ax{2}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, daIFP{di}{bi}{1}{1}, 'Color', 'k', 'LineStyle', ls{1}, 'LineWidth',1.5)
    plot(treal{di}, daIFN{di}{bi}{1}{1}, 'Color', 'k', 'LineStyle', ls{2}, 'LineWidth',1.5)
    fill_between(treal{di}', daIFP{di}{bi}{1}{1}, daIFN{di}{bi}{1}{1} ,  c{bi}, fa(bi));
end
% ylabel("psu m/s", "FontName","Book Antiqua", "Interpreter","latex")
% 

ylabel("m/s", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("$u$: $y$-shear", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = 'sym';
% % opts.XTick = 0;
ax{2} = fix_ax(ax{2}, opts);
ax{2}.XTick = [];



ax{3}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, laIFS{di}{bi}{2}{1}, 'Color', 'k', 'LineStyle', ls{1}, 'LineWidth',1.5)
    plot(treal{di}, laIFB{di}{bi}{2}{1}, 'Color', 'k', 'LineStyle', ls{2}, 'LineWidth',1.5)
    fill_between(treal{di}', laIFS{di}{bi}{2}{1}, laIFB{di}{bi}{2}{1} ,  c{bi}, fa(bi));
end
% ylabel("psu m/s", "FontName","Book Antiqua", "Interpreter","latex")
% 

ylabel("psu", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("$s$: $\sigma$-stratification", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = 'pos';
% % % opts.XTick = 0;
ax{3} = fix_ax(ax{3}, opts);
ax{3}.XTick = [];


ax{4}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, daIFP{di}{bi}{2}{1}, 'Color', 'k', 'LineStyle', ls{1}, 'LineWidth',1.5)
    plot(treal{di}, daIFN{di}{bi}{2}{1}, 'Color', 'k', 'LineStyle', ls{2}, 'LineWidth',1.5)
    fill_between(treal{di}', daIFP{di}{bi}{2}{1}, daIFN{di}{bi}{2}{1} ,  c{bi}, fa(bi));
end
% ylabel("psu m/s", "FontName","Book Antiqua", "Interpreter","latex")
% 

ylabel("psu", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("$s$: $y$-stratification", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = 'pos';
% % % opts.XTick = 0;
ax{4} = fix_ax(ax{4}, opts);
ax{4}.XTick = [];



%Flux
ax{5}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, laIFS{di}{bi}{3}{1}, 'Color', 'k', 'LineStyle', ls{1}, 'LineWidth',1.5)
    plot(treal{di}, laIFB{di}{bi}{3}{1}, 'Color', 'k', 'LineStyle', ls{2}, 'LineWidth',1.5)
    fill_between(treal{di}', laIFS{di}{bi}{3}{1}, laIFB{di}{bi}{3}{1} ,  c{bi}, fa(bi));
end
% ylabel("psu m/s", "FontName","Book Antiqua", "Interpreter","latex")
% 

ylabel("psu m/s", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("$f$: $\sigma$-stratification", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = 'sym';
% % opts.XTick = 0;
ax{5} = fix_ax(ax{5}, opts);
% ax{3}.XTick = [];


ax{6}  = prep_ax();hold on
for bi = bidx
    plot(treal{di}, daIFP{di}{bi}{3}{1}, 'Color', 'k', 'LineStyle', ls{1}, 'LineWidth',1.5)
    plot(treal{di}, daIFN{di}{bi}{3}{1}, 'Color', 'k', 'LineStyle', ls{2}, 'LineWidth',1.5)
    fill_between(treal{di}', daIFP{di}{bi}{3}{1}, daIFN{di}{bi}{3}{1},  c{bi}, fa(bi));
end
% ylabel("psu m/s", "FontName","Book Antiqua", "Interpreter","latex")
% 

ylabel("psu m/s", "FontName", "Book Antiqua", "Interpreter", "latex")
axis tight
title("$f$: $y$-stratification", "FontName","Book Antiqua", "Interpreter", "latex")
opts.ylim = 'sym';
% % opts.XTick = 0;
ax{6} = fix_ax(ax{6}, opts);
% ax{6}.XTick = [];


sgtitle(tit{di}, "FontName","Book Antiqua")

TL.TileSpacing = 'tight';
TL.Padding = 'compact';
%linkaxes([ax{:}],'x')
% TL.XLabel.String = 't [h]';
%     t.YLabel.String = 'z';
% TL.XLabel.Interpreter = 'latex';
TL.YLabel.Interpreter = 'latex';
fontname(TL, "Book Antiqua")
