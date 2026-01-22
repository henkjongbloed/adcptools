function plotMap_xs(DS, scale)

[xr, xl] = utm2ll(mesh.xw, mesh.yw, 31, 'wgs84');

if strcmp(DS, 'NMOMNW15') % ad hoc way of matching adcp / ctd and matlab coordinate systems.
    lat0 = 51.8938250;
    lon0 = 4.3198250;
else
    lat0 = 51.8659640;
    lon0 = 4.332862;
end


hold on
geoplot(xr, xl,'w*', 'MarkerSize',5)
if ~strcmp(scale, 'large')
    geoplot(xr, xl,'w')
end
hold off
end