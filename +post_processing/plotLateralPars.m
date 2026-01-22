function plotLateralPars(U, F)

% Plots geometry
ax = 100; ay = ax;
figure;
npars = U{1}.T.velocity_model.npars  ; % same for all branches
% ntak = numel(U);
idx_0 = [1, npars(1) + 1, sum(npars(1:2)) + 1];
idx_ampl = [(idx_0(1)+1):2:(idx_0(2)-2), (idx_0(2)+1):2:(idx_0(3)-2), (idx_0(3)+1):2:sum(npars)];
idx_phase = idx_ampl + 1;

t = tiledlayout(3,2);



t.XLabel.String = 'y [m]';
t.YLabel.String = 'z [m]';
t.TileSpacing = 'tight';
t.Padding = 'tight';


hold on
for j = 1:numel(U)
    
    nexttile([1,2])

    axis tight
    title([U.BN, ': Subtidal Flow [m/s]'])

    nexttile;

    title('M2 Amplitude [m/s]')
    set(gca,'xticklabel',[])
    axis tight

    nexttile;

    title('M2 Phase [rad]')
    axis tight
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])

    nexttile;

    colorbar;
    title('M4 Amplitude [m/s]')

    axis tight


    nexttile;

    title('M4 Phase [rad]')
    axis tight
    set(gca,'yticklabel',[])
    
    [uxy0, vxy0] = U{j}.xs.sn2xy_pars(F{j}.uc{1, 2}, F{j}.vc{1, 2}); % Subtidal flow
    [uxy1, vxy1] = U{j}.xs.sn2xy_pars(F{j}.uc{2, 2}, F{j}.vc{2, 2}); % Tidal flow
    [uxy2_s, vxy2_s] = U{j}.xs.sn2xy_pars(squeeze(F{j}.uc{3, 2}(:,:,1)), squeeze(F{j}.vc{3, 2}(:,:,1))); % Residual flow t,y,z
    [uxy2_b, vxy2_b] = U{j}.xs.sn2xy_pars(squeeze(F{j}.uc{3, 2}(:,:,end)), squeeze(F{j}.vc{3, 2}(:,:,end))); % Residual flow t,y,z

    u0 = ax*uxy0; v0 = ay*vxy0;
    u1 = ax*(uxy0 + uxy1);  v1 = ay*(vxy0 + vxy1);
    u2_b = ax*(uxy0 + uxy1+ uxy2_b);  v2_b = ay*(vxy0 + vxy1+ vxy2_b);
    u2_s = ax*(uxy0 + uxy1+ uxy2_s);  v2_s = ay*(vxy0 + vxy1+ vxy2_s);


    plot(U{j}.mesh_mean.x_middle(1)  , U{j}.mesh_mean.y_middle(1), 'k*') % Bed (right, positive n)
    plot(U{j}.mesh_mean.x_middle(end)  , U{j}.mesh_mean.y_middle(end), 'k*') % Bed (left, positive n)

    plot(U{j}.mesh_mean.x_middle  , U{j}.mesh_mean.y_middle, 'k') % Transect trajectory mean






    %     uxy = F{j}.uc{1, 3};
    %     vxy = F{j}.vc{1, 3};
    %      = U{j}.xs.sn2xy_pars(F{j}.vc{1, 3});
    %     quiver(U{j}.mesh_mean.x_middle, U{j}.mesh_mean.y_middle, u0, v0, 'k', 'LineWidth', 2) % Subtidal flow
    plot(U{j}.mesh_mean.x_middle + u0, U{j}.mesh_mean.y_middle + v0, 'k', 'LineWidth', 2) % Subtidal flow

    % Tidal Ellipses
    %     plot(U{j}.mesh_mean.x_middle + u1 ,  U{j}.mesh_mean.y_middle + v1, 'b') % Tidal flow, DA
    %     plot(U{j}.mesh_mean.x_middle + u2_b,  U{j}.mesh_mean.y_middle + v2_b, 'r') % Tidal flow, Bottom
    %     plot(U{j}.mesh_mean.x_middle + u2_s,  U{j}.mesh_mean.y_middle + v2_s, 'g') % Tidal flow, Surface
    %

    % Tidal Ellipses: Starting Points
    % Tidal Ellipses
    %     quiver(U{j}.mesh_mean.x_middle + u1(1,:) ,  U{j}.mesh_mean.y_middle + v1(1,:), 'b*') % Tidal flow, DA
    %     quiver(U{j}.mesh_mean.x_middle + u2_b(1,:),  U{j}.mesh_mean.y_middle + v2_b(1,:), 'r*') % Tidal flow, Bottom
    %     quiver(U{j}.mesh_mean.x_middle + u2_s(1,:),  U{j}.mesh_mean.y_middle + v2_s(1,:), 'g*') % Tidal flow, Surface



    %     U{j}.B.plot()
    %     plot3(S{j}.Xp', S{j}.Yp', S{j}.Z' , '.')

end
xlabel('x')
ylabel('y')
title('Tidal Flow')
hold off
end