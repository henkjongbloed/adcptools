function plotPe(dat)

rp = dat.opts.reg_pars;
for i = 1:size(dat.Pe, 1)
    Pe(:,:,i) = reshape(dat.Pe(i,:)',length(rp{3}), length(rp{1}),[]);
end
% regP = dat.regP;
if strcmp(dat.opts.reg_vary, 'coupled')
    [X,Y] = meshgrid(symlog(rp{1}, dat.opts.res_near_zero), symlog(rp{3}, dat.opts.res_near_zero));
end
m = makefigure(24,12);

t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
t = tiledlayout(2,3, TileSpacing="tight", Padding="tight");

t.XLabel.String = '$\lambda_c$';
t.YLabel.String = '$\lambda_d$';
t.XLabel.Interpreter = 'latex';
t.YLabel.Interpreter = 'latex';


tit = {'Total error', 'Relative training error', 'Relative generalization error', '$C_1$: Cell-based continuity (log$_{10}$)', '$C_2$: Continuity across cells (log$_{10}$)', '$C_3$: Parameter coherence (log$_{10}$)', '$C_4$: Gradient consistency (log$_{10}$)', '$A_{\lambda}$: Condition number (log$_{10}$)', 'pcg iterations (log$_{10}$)'};
for fi = 4:9
%     s = subplot(3,3,fi);
    s=nexttile;
    if fi>3 && fi<8
        var = log10(matdivide1(Pe(:,:,fi),Pe(1,1,fi)));
        contourf(X,Y, var,50,'LineStyle','none')
    elseif fi<4
        var = matdivide1(Pe(:,:,fi), Pe(1,1,fi)); 
        hold on
        contourf(X, Y, var, 50, 'LineStyle','none')
        contour(X, Y, var, [1,1], 'LineColor','k')
        hold off
        caxis([1-max(max(abs(var-1))), 1+max(max(abs(var-1)))]);
    else        
        var = log10(Pe(:,:,fi));
        contourf(X, Y, var, 50, 'LineStyle','none')
    end
    colormap(gca, flipud(brewermap(20, 'RdBu')));
    title(tit{fi}, 'interpreter', 'latex')
    idxt = [1, round(size(X,2)/4), round(size(X,2)/2), round(3*size(X,2)/4), round(size(X,2))];
    idyt = [1, round(size(Y,1)/4), round(size(Y,1)/2), round(3*size(Y,1)/4), round(size(Y,1))];
    s.XTick = X(1,idxt);
    s.YTick = Y(idyt,1);
    s.XTickLabel = round(symexp(X(1,idxt), dat.opts.res_near_zero), 1, 'significant');
    s.YTickLabel = round(symexp(Y(idyt,1), dat.opts.res_near_zero), 1, 'significant');
    s.TickLabelInterpreter = 'latex';

    %     xtickformat('%.1f') 

%     if fi==1
%         tlx = symexp(cellfun(@str2num, s.XTickLabel), dat.opts.res_near_zero);
%         tly = symexp(cellfun(@str2num, s.YTickLabel), dat.opts.res_near_zero);
%         tlx_new = arrayfun(@num2str, tlx, 'UniformOutput', 0);
%         tly_new = arrayfun(@num2str, tly, 'UniformOutput', 0);
%     end
%     s.XTickLabel = tlx_new;
%     s.YTickLabel = tly_new;
%     xtickformat('%.1f') 
%     ytickformat('%.1f')%     colormap("jet")
%     s.XTickMode = 
% 
%     xlabel('sB')
%     ylabel('l_c (log)')
%     title(sprintf('||p - p_e||. l_p = %1.1d', Lp(lp_ind(fi))))
    c=colorbar;
    c.TickLabelInterpreter = 'latex';
%     colormap(flipud(hot));     caxis(cap)
end
% sgtitle('Solver performance: Relative to unregularized')
%     subplot(2,4,fi+4)
%     contourf(sB, Lc, log10(squeeze(solB(:,lp_ind(fi), :))),20,'LineWidth',1)
%     xlabel('sB')
%     ylabel('lc (symlog)')
%     title(sprintf('%s%d', 'MSE B estimation. lp:', Lp(lp_ind(fi))))
%     colorbar;
%     colormap(flipud(hot));     caxis(cab)

end