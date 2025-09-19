function ctd = prep_ctd(DS, ctd, tak, V, water_level, xs, mesh_mean, bathy)

if strcmp(DS,'NMOMNW15')    % Change ADCP coordinates rather than CTD coordinates!!
      [ctd{tak}.pos(:,1), ctd{tak}.pos(:,2)]      = wgs2utm(ctd{tak}.Y,ctd{tak}.X, 31, 'N'); %Convert LatLon to UTM %Correct
      time = rescale(ctd{tak}.T, datenum(V.time(1)), datenum(V.time(end)));
      ctd{tak}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
      ctd{tak}.pos(:,3) = water_level.get_water_level(ctd{tak}.t) - ctd{tak}.Z;
%     [ctdH, S.i] = getHour(ctd.T2, hour/2);% + datenum('9-Aug-2014 00:00:00'); % Correct
%     S.t = rescale(ctdH, min(U.adcpH), max(U.adcpH)); % Everything works perfectly
%     % Space
      xnew = center_cloud([ctd{tak}.pos(:,1), ctd{tak}.pos(:,2)], xs.origin'); % Because of offset in CTD coordinates relative to ADCP - quick fix
      ctd{tak}.pos(:,1) = xnew(:,1); ctd{tak}.pos(:,2) = xnew(:,2); % Only needed for 2015 data
      
elseif strcmp(DS,'OMHA14')
    [ctd{tak}.lon, ctd{tak}.lat, ~]      = rd2wgs(ctd{tak}.X, ctd{tak}.Y); %Convert RD to LatLon %Correct
    [ctd{tak}.pos(:,1), ctd{tak}.pos(:,2)]      = wgs2utm(ctd{tak}.lat, ctd{tak}.lon, 31, 'N'); %Convert LatLon to UTM %Correct
    time = rescale(ctd{tak}.T, datenum(V.time(1)), datenum(V.time(end)));
    ctd{tak}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
    ctd{tak}.pos(:,3) = -ctd{tak}.Z;% Upward facing z axis - max(-ctd{tak}.Z - water_level.get_water_level(ctd{tak}.t));
    [ctd{tak}.s, ctd{tak}.n] = xs.xy2sn(ctd{tak}.pos(:,1), ctd{tak}.pos(:,2));
    %     ctddt = 6000;
%     [ctdH, S.i] = getHour(ctd.T, ctddt); % dt = 6000 -> 1 sample per 6s?
     %if strcmp(U.BN, ' Old Meuse North')
%         [S.t, S.i] = manualOMN(ctd.T, U.adcpH);
%     end
end

Xs.t = time;
Xs.x = ctd{tak}.pos(:,1);
Xs.y = ctd{tak}.pos(:,2);
temp.solver.mesh = mesh_mean;
temp.solver.bathy = bathy;
temp.solver.water_level = water_level;
[Hctd, ~, Zbctd]  = get_H(Xs, temp, 0);
hctd = diag(Hctd);
zbctd = diag(Zbctd);
%zpos = ctd{tak}.pos(:,3) + min(ctd{tak}.pos(:,3));
sigctd = (ctd{tak}.pos(:,3)-zbctd)./hctd; % Rescale and scale back to z
sig_corr = rescale(sigctd, 0, 1);
%plotctd(ctd{tak}, mesh_mean, water_level)

ctd{tak}.pos(:,3) = sigctd.*hctd + zbctd;

% figure;
% mesh_mean.plot3()
% hold on
% scatter3(ctd{tak}.pos(:,1), ctd{tak}.pos(:,2), ctd{tak}.pos(:,3))
end