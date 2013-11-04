% SEC2UTC         - convert time in total GPS seconds to UTC time, incl leap sec
%
% Input: time - elapsed seconds since start of GPS time, 1/6/1980
% Output: [yyyy,mm,ss,hh,mm,ss.sss]
%
% See also: GPST2SEC, UTC2GPS, GPS2UTC
% Created: Mike Moreau 4/2002

function [utctime] = sec2utc(time)

utctime = gps2utc(sec2gpst(time));
