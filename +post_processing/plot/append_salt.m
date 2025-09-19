function ax = append_salt(ax, salt, s, xq, yq, ssc, maxs)


s(s<0) = 0;

clr0 = helpers.cmaps('salmap');
clrs = interp1(linspace(0, maxs, 100), clr0, 0:1:maxs);

% prepare for area plots
[sx, sy] = salt.solver.xs.sn2xy_vel(s, zeros(size(s)));

sxp = xq +  ssc*sx; sxm = xq -  ssc*sx;
syp = yq +  ssc*sy; sym = yq -  ssc*sy;
% disp(clrs(round(mean(s)) + 1,:))
% disp(["Avg salinity: ", num2str(mean(s))])
fill(ax,[sxp, fliplr(sxm)], [syp, fliplr(sym)], clrs(round(max(s)) + 1,:))

end