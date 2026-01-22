function plotTe(dat)

rp = dat.opts.reg_pars;
for i = 1:size(dat.Pe, 1)
    Pe(:,:,i) = reshape(dat.Pe(i,:)',length(rp{3}), length(rp{1}),[]);
end
% regP = dat.regP;
if strcmp(dat.opts.reg_vary, 'coupled')
    [X,Y] = meshgrid(symlog(rp{1}, dat.opts.res_near_zero), symlog(rp{3}, dat.opts.res_near_zero));
end


m=makefigure(8,5);
% tit = {'Total error', 'Training Error', 'Test error', 'Continuity: Inside cell (log)', 'Continuity: Inter-cell (log)', 'Coherence (log)', 'Consistency (log)'};
% for fi = 1:length(tit)
s = subplot(1,1,1);
fi=3;
if fi>3
    var = log10(matdivide1(Pe(:,:,3),Pe(1,1,fi)));
else
    var = matdivide1(Pe(:,:,fi),Pe(1,1,fi));
end

[minvar, I] = min(var(:));
[mi, mj] = ind2sub(size(var), I);
hold on
contourf(X,Y, var,50,'LineStyle','none')
contour(X,Y, var, [1,1], 'LineColor','k')
colormap(gca, flipud(brewermap(20, 'RdBu')));
caxis([1-max(max(abs(var-1))), 1+max(max(abs(var-1)))]);
hold off
tlx = symexp(cellfun(@str2num, s.XTickLabel), dat.opts.res_near_zero);
tly = symexp(cellfun(@str2num, s.YTickLabel), dat.opts.res_near_zero);
tlx_new = arrayfun(@num2str, tlx, 'UniformOutput', 0);
tly_new = arrayfun(@num2str, tly, 'UniformOutput', 0);

s.XTickLabel = tlx_new;
s.YTickLabel = tly_new;
xt = get(gca, 'XAxis');
xt.TickLabelInterpreter = 'latex'; % latex for x-axis
yt = get(gca, 'YAxis');
yt.TickLabelInterpreter = 'latex'; % latex for x-axis
c=colorbar;
c.TickLabelInterpreter = 'latex';

xlabel('$\lambda_c$', 'interpreter', 'latex')
ylabel('$\lambda_d$', 'interpreter', 'latex')
title('Relative validation error', 'interpreter', 'latex')
%     colormap("jet")

%
%     xlabel('sB')
%     ylabel('l_c (log)')
%     title(sprintf('||p - p_e||. l_p = %1.1d', Lp(lp_ind(fi))))
%     colormap(flipud(hot));     caxis(cap)
% end
% sgtitle('Solver performance: Relative to unregularized')
%     subplot(2,4,fi+4)
%     contourf(sB, Lc, log10(squeeze(solB(:,lp_ind(fi), :))),20,'LineWidth',1)
%     xlabel('sB')
%     ylabel('lc (symlog)')
%     title(sprintf('%s%d', 'MSE B estimation. lp:', Lp(lp_ind(fi))))
%     colorbar;
%     colormap(flipud(hot));     caxis(cab)

end