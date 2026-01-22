function plotMap_ctd(sal)

[xr, xl] = utm2ll(sal.solver.mesh.xw, sal.solver.mesh.yw, 31, 'wgs84');


hold on
geoplot(xr, xl,'w*', 'MarkerSize',5)
if ~strcmp(scale, 'large')
    geoplot(xr, xl,'w')
end
hold off
end