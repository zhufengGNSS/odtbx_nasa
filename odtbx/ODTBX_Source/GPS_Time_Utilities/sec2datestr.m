% SEC2DATESTR      - convert time in total GPS seconds to a formated date and time string
%
% [utcstr] = sec2datestr(time,output_format)
% 
% Input:    time - elapsed seconds since start of GPS time, 1/6/1980
%           output_format: 0-   6/1/2002 12:00:00  (with zero padding on times)
%                          1-   1-JUN-2002 12:00:00  (with zero padding on times)
%                          2-   1 Jun 2002 12:00:00.000000000  (with zero padding on times)
%                          3-   YYYYMMDD.HHMMSS  (with zero padding)
% Output:   datestr:       formatted string
% 
% See also: UTC2DATESTRING,CALENDAR2DATASTR
% Written: Mike Moreau 4/2002

function [utcstr] = sec2datestr(time,output_format)

if nargin <2
    output_format = 0;
elseif (output_format ~= 0) & (output_format ~= 1) & (output_format ~= 2) & (output_format ~= 3)
    error('Invalid output format.')
end

utctime = gps2utc(sec2gpst(time(1)));

utcstr = calendar2datestr(utctime,output_format);