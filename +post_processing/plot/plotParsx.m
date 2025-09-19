function plotParsx(U)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
figure;
% npars = U.T.velocity_model.npars  ;
names = U.T.velocity_model.names;
% ntak = numel(U);
idx_0 = [find(strcmp(names, 'u0: M0A')) , ...
        find(strcmp(names, 'v0: M0A')) , ...
        find(strcmp(names, 'w0: M0A'))];

idx_a2 = [find(strcmp(names, 'u0: M2A'))];
idx_b2 = [find(strcmp(names, 'u0: M2B'))];

idx_a4 = [find(strcmp(names, 'u0: M4A'))];
idx_b4 = [find(strcmp(names, 'u0: M4B'))];


t = tiledlayout(3,2);

nexttile([1,2])
U.mesh_mean.plot(U.pars(:,idx_0))
amax = max(max(abs(U.pars(:,idx_0)), [], 'omitnan'));
colormap(gca, U.velmap)
caxis([-amax, amax])
colorbar;
axis tight
title([U.BN, ': Subtidal Flow [m/s]'])
set(gca, 'XDir','reverse') % Very important

nexttile;
tempa2 = sqrt(U.pars(:,idx_a2).^2 + U.pars(:,idx_b2).^2);
tempp2 = atan(U.pars(:,idx_b2)./U.pars(:,idx_a2));

tempa4 = sqrt(U.pars(:,idx_a4).^2 + U.pars(:,idx_b4).^2);
tempp4 = atan(U.pars(:,idx_b4)./U.pars(:,idx_a4));
U.mesh_mean.plot(tempa2);
 
amax = max(max(abs(tempa2), [], 'omitnan'));
colormap(gca, U.velmap)
caxis([0, amax])
colorbar;
title('M2 Amplitude [m/s]')
set(gca,'xticklabel',[])
axis tight
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(tempp2)
% amax = max(max(abs(U.pars(:,idx_phase)), [], 'omitnan'));
colormap(gca, U.phimap)
caxis([-pi/2, pi/2])
colorbar;
title('M2 Phase [rad]')
axis tight
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca, 'XDir','reverse') % Very important

nexttile;
U.mesh_mean.plot(tempa4)
amax = max(max(abs(tempa4), [], 'omitnan'));
colormap(gca, U.velmap)
caxis([0, amax])
colorbar;
title('M4 Amplitude [m/s]')
set(gca, 'XDir','reverse') % Very important
axis tight


nexttile;
U.mesh_mean.plot(tempp4)
% amax = max(max(abs(U.pars(:,idx_phase)), [], 'omitnan'));
colormap(gca, U.phimap)
% colormap(gca, U.velmap)
caxis([-pi/2, pi/2])
colorbar;
title('M4 Phase [rad]')
axis tight
set(gca,'yticklabel',[])
set(gca, 'XDir','reverse') % Very important

t.XLabel.String = 'y [m]';
t.YLabel.String = 'z [m]';
t.TileSpacing = 'tight';
t.Padding = 'tight';

end