function [tc, idxtc] = clip_time(t, tlim, dt)
% input: datetime vectors and limits

% output: clipped datetime array and original indices
% dt: time in days that the limits need to be enlarged with (symmetric)

tn = datenum(t);
tlimn = datenum(tlim);
tlimn(1) = tlimn(1) - dt;
tlimn(2) = tlimn(2) + dt;

[tcn, idxtc] = clip(tn, tlimn);

tc = datetime(tcn, 'ConvertFrom', 'datenum');

end