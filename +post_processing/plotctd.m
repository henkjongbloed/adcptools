function plotctd(struc, mesh, wl)

% t=1;
% if strcmp(DS,'NMOMNW15')    % Change ADCP coordinates rather than CTD coordinates!!
%       [ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2)]      = wgs2utm(ctd{tak(t)}.Y,ctd{tak(t)}.X, 31, 'N'); %Convert LatLon to UTM %Correct
%       time = rescale(ctd{tak(t)}.T, datenum(V.time(1)), datenum(V.time(end)));
%       ctd{tak(t)}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
%       ctd{tak(t)}.pos(:,3) = water_level.get_water_level(ctd{tak(t)}.t) - ctd{tak(t)}.Z;
% %     [ctdH, S.i] = getHour(ctd.T2, hour/2);% + datenum('9-Aug-2014 00:00:00'); % Correct
% %     S.t = rescale(ctdH, min(U.adcpH), max(U.adcpH)); % Everything works perfectly
% %     % Space
% %     xnew = center_cloud([S.Xold, S.Yold], U.xs.origin'); % Because of offset in CTD coordinates relative to ADCP - quick fix
% %     S.X = xnew(:,1); S.Y = xnew(:,2); % Only needed for 2015 data
%       
% elseif strcmp(DS,'OMHA14')
%     [ctd{tak(t)}.lon, ctd{tak(t)}.lat, ~]      = rd2wgs(ctd{tak(t)}.X, ctd{tak(t)}.Y); %Convert RD to LatLon %Correct
%     [ctd{tak(t)}.pos(:,1), ctd{tak(t)}.pos(:,2)]      = wgs2utm(ctd{tak(t)}.lat, ctd{tak(t)}.lon, 31, 'N'); %Convert LatLon to UTM %Correct
%     time = rescale(ctd{tak(t)}.T, datenum(V.time(1)), datenum(V.time(end)));
%     ctd{tak(t)}.t = datetime(time, 'ConvertFrom', 'datenum'); % Rescale time such that measurements fall within experiment time
%     ctd{tak(t)}.pos(:,3) = water_level.get_water_level(ctd{tak(t)}.t) - (ctd{tak(t)}.Z - min(ctd{tak(t)}.Z));
% %     ctddt = 6000;
% %     [ctdH, S.i] = getHour(ctd.T, ctddt); % dt = 6000 -> 1 sample per 6s?
%      %if strcmp(U.BN, ' Old Meuse North')
% %         [S.t, S.i] = manualOMN(ctd.T, U.adcpH);
% %     end
% end


figure;
subplot(4,1,1)
plot(struc.t, struc.pos(:,3), "*")
subplot(4,1,2)
plot(struc.t, struc.Z, "*")
subplot(4,1,3)
plot(struc.t, wl.get_water_level(struc.t), "*")
subplot(4,1,4)
hold on
%mz = 
scatter(struc.t, -struc.Z - max(-struc.Z - wl.get_water_level(struc.t)), 10*struc.S)
plot(struc.t, wl.get_water_level(struc.t))