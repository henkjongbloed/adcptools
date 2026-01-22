function ax = plot_dynamics_residual(di, bidx, flow, salt,D,  AF, AIFB, AIFS)
%% Plot subtidal dynamics with top view

c = get_c();
%

maxs{1} = 3; % max depth-averaged salinities -> integers
maxs{2} = 18;
% figure;

% Preliminary plot
rp = 1;
% One figure per dataset, one subplot per time
% close all
yc = 2:4:(numel(D{di}{bidx(1)}.zb)-1); % coarsening of lateral profiles.
yca = 1:2:numel(D{di}{bidx(1)}.zb  );
ax=nexttile;
ti=1;

ssc = 10;
if di ==1
    fsc = 50*ssc; % flow scale
else
    fsc = 100*ssc; % flow scale
end

for bi = bidx

    xb_ = [flow{di}{bi}.solver.mesh.x_middle(1), flow{di}{bi}.solver.mesh.x_middle(end)];
    yb_ = [flow{di}{bi}.solver.mesh.y_middle(1), flow{di}{bi}.solver.mesh.y_middle(end)];

    xq = linspace(xb_(1), xb_(2), D{di}{bi}.sz(2));
    yq = linspace(yb_(1), yb_(2), D{di}{bi}.sz(2));

    hold on
    var = 2; %salt
    s0 = AF{di}{bi}{var}{rp,1}{3}(ti, :);
    sb = AIFB{di}{bi}{var}{rp,1}(ti, :);
    ss = AIFS{di}{bi}{var}{rp,1}(ti, :);


    
    ax = append_salt(ax, salt{di}{bi}, sb, xq, yq, ssc, maxs{di});
    ax = append_salt(ax, salt{di}{bi}, s0, xq, yq, ssc, maxs{di});
    ax = append_salt(ax, salt{di}{bi}, ss, xq, yq, ssc, maxs{di});
    %             colorbar;


    var = 1; % flow
    % subtidal: index 3 in AF
    [u0, v0] = flow{di}{bi}.solver.xs.sn2xy_vel(AF{di}{bi}{var}{rp,1}{3}(ti, :), AF{di}{bi}{var}{rp, 2}{3}(ti, :)); % depth-averaged subtidal
    [ub, vb] = flow{di}{bi}.solver.xs.sn2xy_vel(AIFB{di}{bi}{var}{rp,1}(ti, :), AIFB{di}{bi}{var}{rp,2}(ti, :));
    [us, vs] = flow{di}{bi}.solver.xs.sn2xy_vel(AIFS{di}{bi}{var}{rp,1}(ti, :), AIFS{di}{bi}{var}{rp,2}(ti, :));
    % flow
    if di ==1
        quiver(ax, xq(yc+2), yq(yc+2), fsc*us(yc+2), fsc*vs(yc+2), 'b', 'LineWidth', 2, 'AutoScale', 'off', 'Color', 'b') % Subtidal flow surf
        quiver(ax, xq(yca), yq(yca), fsc*u0(yca), fsc*v0(yca), 'k', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
        quiver(ax, xq(yc), yq(yc), fsc*ub(yc), fsc*vb(yc), 'r', 'LineWidth', 2, 'AutoScale', 'off', 'Color', 'r') % Subtidal flow bottom
    else
        quiver(ax, xq(yc), yq(yc), fsc*ub(yc), fsc*vb(yc), 'r', 'LineWidth', 2, 'AutoScale', 'off', 'Color', 'r') % Subtidal flow bottom
        quiver(ax, xq(yc+2), yq(yc+2), fsc*us(yc+2), fsc*vs(yc+2), 'b', 'LineWidth', 2, 'AutoScale', 'off', 'Color', 'b') % Subtidal flow surf
        quiver(ax, xq(yca), yq(yca), fsc*u0(yca), fsc*v0(yca), 'k', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
    end

%     arrowSurfLength = sqrt(us(yc+4).^2+vs(yc+4).^2)% for visual inspection of quiver arrow scales.

    saltBotLength = max(sb)% for visual inspection of quiver arrow scales.

    hold on
    plot(ax, xb_, yb_, "*k")
    plot(ax, xq, yq, "k")
    hold off
end
axis square

end
