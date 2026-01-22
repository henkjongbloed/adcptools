function [ax, maxuu, maxss] = plot_dynamics_instantaneous(ax, flow, salt, D, AF, IFB, IFS, di, bi, ti)

% Preliminary plot
rp = 1;
% One figure per dataset, one subplot per time
% close all
% close all
yc = 2:9:(numel(D{di}{bi}.zb)-1); % coarsening of lateral profiles.
yca = 1:3:numel(D{di}{bi}.zb  );
maxs{1} = 6; % max depth-averaged salinities -> integers
maxs{2} = 23;

ssc = 10;
if di ==1
    fsc = 20*ssc; % flow scale
else
    fsc = 30*ssc; % flow scale
end

% ax = nexttile;

xb_ = [flow{di}{bi}.solver.mesh.x_middle(1), flow{di}{bi}.solver.mesh.x_middle(end)];
yb_ = [flow{di}{bi}.solver.mesh.y_middle(1), flow{di}{bi}.solver.mesh.y_middle(end)];

xq = linspace(xb_(1), xb_(2), D{di}{bi}.sz(2));
yq = linspace(yb_(1), yb_(2), D{di}{bi}.sz(2));

hold on
var = 2; %salt
s0 = AF{di}{bi}{var}{rp,1}{7}(ti, :);
sb = IFB{di}{bi}{var}{rp,1}(ti, :);
ss = IFS{di}{bi}{var}{rp,1}(ti, :);


% ssc = 20;
ax = append_salt(ax, salt{di}{bi}, sb, xq, yq, ssc, maxs{di});
ax = append_salt(ax, salt{di}{bi}, s0, xq, yq, ssc, maxs{di});
ax = append_salt(ax, salt{di}{bi}, ss, xq, yq, ssc, maxs{di});



var = 1; % flow
[u0, v0] = flow{di}{bi}.solver.xs.sn2xy_vel(AF{di}{bi}{var}{rp,1}{7}(ti, :), AF{di}{bi}{var}{rp, 2}{7}(ti, :)); % depth-averaged tidal
[ub, vb] = flow{di}{bi}.solver.xs.sn2xy_vel(IFB{di}{bi}{var}{rp,1}(ti, :), IFB{di}{bi}{var}{rp, 2}(ti, :));
[us, vs] = flow{di}{bi}.solver.xs.sn2xy_vel(IFS{di}{bi}{var}{rp,1}(ti, :), IFS{di}{bi}{var}{rp, 2}(ti, :));
% flow
% fsc = 30*ssc;
quiver(ax, xq(yca), yq(yca), fsc*u0(yca), fsc*v0(yca), 'k', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
quiver(ax, xq(yc+3), yq(yc+3), fsc*us(yc+3), fsc*vs(yc+3), 'b', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow
quiver(ax, xq(yc), yq(yc), fsc*ub(yc), fsc*vb(yc), 'r', 'LineWidth', 2, 'AutoScale', 'off') % Subtidal flow

maxuu = max(sqrt(us.^2+vs.^2));% for visual inspection of quiver arrow scales.
maxss = max(sb);% for visual inspection of quiver arrow scales.

hold on
plot(ax, xb_, yb_, "*k")
plot(ax, xq, yq, "k")
hold off
end