function plotSubtidalSalt(tak, F, U)
figure;

for j = 1:length(tak)
    s = F{j}.s{1};
    for i = 2:length(F{j}.s) %very ugly
        s = s + F{j}.s{i};
    end
    s = s/length(F{j}.s);
    %npars = U{j}.T.velocity_model.npars  ;
    %subtidal_idx = [1, npars(1) + 1, sum(npars(1:2)) + 1];
    %subtidal_u = U{j}.tid_pars(:,subtidal_idx);
    %     b = tak(j);
    subplot(length(tak), 1, j)
    %lnorm = sqrt(subtidal_u(:,2).^2 + subtidal_u(:,3).^2);
    U{j}.mesh_mean.plot(s)
%     max(subtidal_u(:,2))
%     max(subtidal_u(:,3))

%     hold on
%     quiver(.9*max(U{j}.mesh_mean.n_right), min(U{j}.mesh_mean.z_bottom_right), 1, 0, 1,'k')
%     quiver(.9*max(U{j}.mesh_mean.n_right), min(U{j}.mesh_mean.z_bottom_right), 0, 1, 1,'k')
%     hold off
    colormap(gca, U{j}.salmap)
    c=colorbar;
    
    caxis([0, max(s, [], 'all')])
    c.Label.String = 'Salinity [psu]';
    axis tight
    title(strcat('Salinity: ', U{j}.BN), 'FontSize',12)
    ylabel('z [m]')
    if j==length(tak)
        xlabel('y [m]')
    end
    set(gca, 'XDir','reverse')
end
% print('salt', '-dpng')
end