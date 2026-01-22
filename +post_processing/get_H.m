function [H, Wl, Zb] = get_H(X, sol, xs_proj)


if xs_proj
    [~, yvec] = sol.solver.xs.sn2xy(zeros(size(X.y)), X.y); % convert to x,y -> TODO Enter correct s,n coord!!
else
%     xvec = X.x;
    yvec = X.y; % convert to x,y -> TODO Enter correct s,n coord!!
end

% zvec = sol.solver.bathy.get_bed_elev(xvec, yvec);  % get (time-indep) bathy %DEPRECATED
zb = n2zb(yvec, sol.solver.mesh);
wl = sol.solver.bathy.water_level.get_water_level(datetime(X.t, 'ConvertFrom', 'datenum'));

[Wl, Zb] = ndgrid(wl, zb);

H = Wl - Zb;

end