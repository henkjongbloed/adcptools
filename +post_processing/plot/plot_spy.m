function plot_spy = plot_spy(dat)


makefigure(15,10)
t = tiledlayout(2,3);
A = dat.M'*dat.M + dat.C'*dat.C + dat.Cg'*dat.Cg + ...
    dat.D'*dat.D + dat.D2'*dat.D2;
mat = {dat.M'*dat.M, dat.C'*dat.C, dat.Cg'*dat.Cg,...
    dat.D'*dat.D, dat.D2'*dat.D2, A};
tit = {'$M^TM$', '$C_1^TC_1$', '$C_2^TC_2$'...
    '$C_3^TC_3$', '$C_4^TC_4$', '$A$'};
for i=1:length(tit)
nexttile;
spy(mat{i});
title(tit{i}, 'interpreter', 'latex')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
end
end