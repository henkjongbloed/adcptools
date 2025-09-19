function plotNoiseLD(dat, lc_idx, est_opts)

rp = dat.opts.reg_pars;
SA = reshape(dat.avg_err, length(rp{3}),length(rp{1}),[]);
if strcmp(dat.opts.reg_vary, 'coupled')
    [X,Y] = meshgrid(est_opts.noise_levels, symlog(rp{3}, dat.opts.res_near_zero));
end
% figure('units','normalized','outerposition',[0 0 1 1]);
% tit = {'Total error', 'Training Error', 'Test error', 'Continuity: Inside cell', 'Continuity: Inter-cell', 'Coherence', 'Consistency'};
t = tiledlayout('flow', TileSpacing="tight", Padding="tight");

for fi = lc_idx
    s=nexttile;
    var = squeeze(SA(:,fi,:));
    contourf(X',Y', log10(var'), 50, 'LineStyle','none')
%     xlabel('$\sigma$', 'interpreter', 'latex')
% ylabel('$\lambda_d$', 'interpreter', 'latex')
    %     if fi==1
    %         tlx = symexp(cellfun(@str2num, s.XTickLabel), dat.opts.res_near_zero);
    %         tly = symexp(cellfun(@str2num, s.YTickLabel), dat.opts.res_near_zero);
    %         tlx_new = arrayfun(@num2str,tlx, 'UniformOutput', 0);
    %         tly_new = arrayfun(@num2str, tly, 'UniformOutput', 0);
    %     end
    %     s.XTickLabel = tlx_new;
    %     s.YTickLabel = tly_new;

%     idxt = [1, max([2,round(size(X,2)/4)]), round(size(X,2)/2), round(3*size(X,2)/4), round(size(X,2))];
    idyt = [1, max([2,round(size(Y,1)/4)]), round(size(Y,1)/2), round(3*size(Y,1)/4), round(size(Y,1))];
%     s.XTick = X(1,idxt);
    s.YTick = Y(idyt,1);
%     s.XTickLabel = round(symexp(X(1,idxt), dat.opts.res_near_zero), 1, 'significant');
    s.YTickLabel = round(symexp(Y(idyt,1), dat.opts.res_near_zero), 1, 'significant');
    s.TickLabelInterpreter = 'latex';

    title(sprintf('log$_{10}$ error: $\\lambda_c$ = %.2f', rp{1}(fi)), 'interpreter', 'latex')
%     disp(rp{1}(fi))
    %     colormap("jet") 
    colormap(gca, flipud(brewermap(20, 'RdBu')));
    c=colorbar;   
    c.TickLabelInterpreter = 'latex';
%     colorbar;
end
t.XLabel.String = '$\sigma$';
t.YLabel.String = '$\lambda_d$';
t.XLabel.Interpreter = 'latex';
t.YLabel.Interpreter = 'latex';
% sgtitle('Sensitivity analysis: Additive noise')

% for ni = noise_ind
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     SAP = dat.avg_perr;
%     for fi = 1:size(SAP)
%         s=subplot(6,10,fi);
%         sap = reshape(squeeze(SAP(fi,:,:)),length(rp{3}), length(rp{1}),[]);
%         contourf(X,Y,matdivide1(sap(:,:,ni),sap(1,1,ni)),20,'LineWidth',1)
%         %     xlabel('lc')
%         %     ylabel('ld')
%         title(names{fi})
%         if fi==1
%             tlx = symexp(cellfun(@str2num, s.XTickLabel), dat.opts.res_near_zero);
%             tly = symexp(cellfun(@str2num, s.YTickLabel), dat.opts.res_near_zero);
%             tlx_new = arrayfun(@num2str, tlx, 'UniformOutput', 0);
%             tly_new = arrayfun(@num2str, tly, 'UniformOutput', 0);
%         end
%         s.XTickLabel = tlx_new;
%         s.YTickLabel = tly_new;
%         colorbar;
%     end
%     sgtitle(sprintf('Sensitivity analysis per parameter: noise level $\\sigma$ = %i', est_opts.noise_levels(ni)), 'interpreter', 'latex')
% end

end