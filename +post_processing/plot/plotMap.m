function plotMap(DS, scale)

if strcmp(DS, 'NMOMNW15')
    lat0 = 51.8938250;
    lon0 = 4.3198250;
else
    lat0 = 51.8659640;
    lon0 = 4.332862;
end

if strcmp(scale, 'large')
    latlim = lat0 + 6*.006*[-.5, 4];
    lonlim = lon0 + 10*.0025*[-15, 1];
else
    latlim = lat0 + .006*[-1, 1];
    lonlim = lon0 + + .0025*[-1, 1];
end

geolimits(latlim ,lonlim)
geobasemap("satellite")

end