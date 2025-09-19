function [H, Wl, Zvec] = get_H_mesh(X, sol)


xvec = zeros(size(X.y));%X.x;
yvec = X.y; % convert to x,y -> TODO Enter correct s,n coord!!


%zvec = sol.solver.bathy.get_bed_elev(xvec, yvec);  % get (time-indep) bathy
zvec = n2zb(yvec, sol.solver.mesh);
wl = sol.solver.bathy.water_level.get_water_level_model(datetime(X.t, 'ConvertFrom', 'datenum'));

[Wl, Zvec] = ndgrid(wl, zvec);

H = Wl - Zvec;

end