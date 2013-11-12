% UTC2SEC         - convert calendar date in UTC to time in total elapsed seconds since GPS epoch 1980
%
% Input: [yyyy,mm,ss,hh,mm,ss.sss]
% Output: time - elapsed seconds since start of GPS time, 1/6/1980
%
% See also: GPST2SEC, UTC2GPS, GPS2UTC
% Created: Mike Moreau 4/2002

function [time] = utc2sec(utctime)

time = gpst2sec(utc2gps(utctime));
