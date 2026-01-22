function plotSubtidalFlow(tak, U)
figure;


for j = 1:length(tak)
    npars = U{j}.T.velocity_model.npars  ;
    idx_0 = [1, npars(1) + 1, sum(npars(1:2)) + 1];
    subtidal_u = U{j}.tid_pars(:,idx_0);
    %     b = tak(j);
    subplot(length(tak), 1, j)
%     lnorm = 1; %sqrt(subtidal_u(:,2).^2 + subtidal_u(:,3).^2);
    U{j}.mesh_mean.plot(subtidal_u)
%     max(subtidal_u(:,2))
%     max(subtidal_u(:,3))

    hold on
%     quiver(.9*max(U{j}.mesh_mean.n_right), min(U{j}.mesh_mean.z_bottom_right), 1, 0, 1,'k')
%     quiver(.9*max(U{j}.mesh_mean.n_right), min(U{j}.mesh_mean.z_bottom_right), 0, 1, 1,'k')
    hold off
    colormap(gca, U{j}.velmap)
    caxis([-.7, .7])
    c=colorbar;
    c.Label.String = 'Flow [m/s]';

    axis tight
    title(strcat('Flow: ', U{j}.BN), 'FontSize',12)
    ylabel('z [m]')
    if j==length(tak)
        xlabel('y [m]')
    end
    set(gca, 'XDir','reverse') % Very important
end

% print('flow', '-dpng')
end