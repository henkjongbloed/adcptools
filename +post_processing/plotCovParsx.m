function plotCovParsx(U)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
figure;
npars = U.T.velocity_model.npars  ;
% ntak = numel(U);
idx_0 = [1, npars(1) + 1, sum(npars(1:2)) + 1];
idx_ampl = [(idx_0(1)+1):2:(idx_0(2)-2), (idx_0(2)+1):2:(idx_0(3)-2), (idx_0(3)+1):2:sum(npars)];
idx_phase = idx_ampl + 1;

t = tiledlayout(3,2);

for i = 1:size(U.cov_tid_pars,1) %Prepare plotted data to check covariances
    c(i,1) = trace(squeeze(U.cov_tid_pars(i,:,:)));
    c(i,2) = U.cov_tid_pars(i,idx_0(1),idx_0(1));
    c(i,3) = U.cov_tid_pars(i,idx_0(2),idx_0(2));
    c(i,4) = U.cov_tid_pars(i,idx_0(3),idx_0(3));
% M2
    c(i,5) = U.cov_tid_pars(i,idx_ampl(1),idx_ampl(1));
    c(i,6) = U.cov_tid_pars(i,idx_ampl(3),idx_ampl(3));
    c(i,7) = U.cov_tid_pars(i,idx_ampl(5),idx_ampl(5));

    c(i,8) = U.cov_tid_pars(i,idx_phase(1),idx_phase(1));
    c(i,9) = U.cov_tid_pars(i,idx_phase(3),idx_phase(3));
    c(i,10) = U.cov_tid_pars(i,idx_phase(5),idx_phase(5));
% M4
    c(i,11) = U.cov_tid_pars(i,idx_ampl(2),idx_ampl(2));
    c(i,12) = U.cov_tid_pars(i,idx_ampl(4),idx_ampl(4));
    c(i,13) = U.cov_tid_pars(i,idx_ampl(6),idx_ampl(6));

    c(i,14) = U.cov_tid_pars(i,idx_phase(2),idx_phase(2));
    c(i,15) = U.cov_tid_pars(i,idx_phase(4),idx_phase(4));
    c(i,16) = U.cov_tid_pars(i,idx_phase(6),idx_phase(6));
end

nexttile;

U.mesh_mean.plot(c(:,1))
% amax = max(max(abs(U.tid_pars(:,idx_0)), [], 'omitnan'));
colormap(gca, U.velmap)
%caxis([-amax, amax])
colorbar;
axis tight
title([U.BN, ': Trace of Covariance Matrix'])
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(c(:,2:4));
% amax = max(max(abs(U.tid_pars(:,idx_ampl(:,1))), [], 'omitnan'));
colormap(gca, U.velmap)
%caxis([0, amax])
colorbar;
title('Subtidal variances')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(c(:,5:7));
% amax = max(max(abs(U.tid_pars(:,idx_ampl(:,1))), [], 'omitnan'));
colormap(gca, U.velmap)
%caxis([0, amax])
colorbar;
title('M2 Ampl variances')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(c(:,8:10));
% amax = max(max(abs(U.tid_pars(:,idx_ampl(:,1))), [], 'omitnan'));
colormap(gca, U.velmap)
%caxis([0, amax])
colorbar;
title('M2 Phase variances')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(c(:,11:13));
% amax = max(max(abs(U.tid_pars(:,idx_ampl(:,1))), [], 'omitnan'));
colormap(gca, U.velmap)
%caxis([0, amax])
colorbar;
title('M4 Ampl variances')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(c(:,14:16));
% amax = max(max(abs(U.tid_pars(:,idx_ampl(:,1))), [], 'omitnan'));
colormap(gca, U.velmap)
%caxis([0, amax])
colorbar;
title('M4 Phase variances')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important


t.XLabel.String = 'y [m]';
t.YLabel.String = 'z [m]';
t.TileSpacing = 'tight';
t.Padding = 'tight';

end