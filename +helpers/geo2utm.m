function [UTMx, UTMy, zone] = geo2utm(lat,long,zone)
%
%    Copyright 2009,2010 Bart Vermeulen, Maximiliano Sassi
%
%    This file is part of ADCPTools.
%
%    ADCPTools is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    ADCPTools is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ADCPTools.  If not, see <http://www.gnu.org/licenses/>.

%% find utm zone for given data
% select good data
goodff=(~isnan(lat)) & (~isnan(long));
if ~any(goodff)
    UTMx=long;
    UTMy=lat;
    zone='';
    return; 
end

if nargin() < 3 || isempty(zone)
    
    % find zone
    mlat=mean(lat,'omitnan');
    mlong=mean(long,'omitnan');
    lts = [-80:8:72 84]';
    lns = (-180:6:180)';
    latzones = char([67:72 74:78 80:88]');
    
    indx = find(lts <= mlat);
    ltindx = indx(max(indx));
    
    indx = find(lns <= mlong);
    lnindx = indx(max(indx));
    
    if ltindx < 1 || ltindx > 21
        ltindx = [];
    elseif ltindx == 21
        ltindx = 20;
    end
    
    if lnindx < 1 || lnindx > 61
        lnindx = [];
    elseif lnindx == 61
        lnindx = 60;
    end
    
    zone = [num2str(lnindx) latzones(ltindx)];
    
end % if isempty(zone)

% transform with suitable package
% if exist('utm_fwd','file')==2 % GeographicLib is available
%     disp('Using GeographicLib for UTM forward transformation')
%     if mean(lat(goodff))>0, northp=true; else northp=false; end
%     [UTMx,UTMy]=utm_fwd(lnindx,northp,lat,long);
% elseif license('checkout','map_toolbox')
%     disp('Using Mapping Toolbox for UTM forward transformation')
% %     zone=utmzone(lat(goodff),long(goodff));
%     ellipsoid = utmgeoid(zone);
%     adcpmap = defaultm('utm');
%     adcpmap.zone = zone;
%     adcpmap.geoid = ellipsoid;
%     adcpmap.flatlimit = [];
%     adcpmap.maplatlimit = [];
%     adcpmap = defaultm(adcpmap);
%     [UTMx,UTMy]=mfwdtran(adcpmap,lat,long);
% else
%     disp('Using custom UTM forward transformation')

size_in=size(lat);

lns = (-180:6:180)';
indx = find(lns <= mean(long(goodff)));
zone = indx(max(indx));
if zone < 1 || zone > 61
    zone = [];
elseif zone == 61
    zone = 60;
end
long0=(zone*6-183)/180*pi;
long=reshape(long/180*pi,[],1);
lat=reshape(lat/180*pi,[],1);
a=6378.137;
finv=298.257223563; f=1/finv;
if mean(lat(goodff))>0
    N0=0;
else
    N0=10000; 
end
k0=0.9996;
E0=500;
n=f/(2-f);
t=sinh(atanh(sin(lat))-2*sqrt(n)/(1+n)*atanh(2*sqrt(n)*sin(lat)/(1+n)));
zeta=atan(t./cos(long-long0));
eta=atanh(sin(long-long0)./sqrt(1+t.^2));
A=a/(1+n)*(1+n.^2/4+n^4/64);
alph(1)=1/2*n-2/3*n^2+5/16*n^3;
alph(2)=13/48*n^2-3/5*n^3;
alph(3)=61/240*n^3;
jvec=1:3;
UTMx=E0+k0*A*(eta+sum(bsxfun(@times,alph,cos(2*bsxfun(@times,jvec,zeta)).*sinh(2*bsxfun(@times,jvec,eta))),2));
UTMy=N0+k0*A*(zeta+sum(bsxfun(@times,alph,sin(2*bsxfun(@times,jvec,zeta)).*cosh(2*bsxfun(@times,jvec,eta))),2));
UTMx=UTMx*1000;
UTMy=UTMy*1000;
UTMx=reshape(UTMx, size_in);
UTMy=reshape(UTMy, size_in);
% end

end % function geo2utm

