% GPS2MJD         - Convert GPS time to Modified Julian Date
%
%   [mjd,sod] = gps2mjd(week,sec)
%
%  Converts GPS TIME in the format:
%     [week,seconds of week] [n,2]
%  to a UTC TIME in modified Julian Date format:
%     [MJD,seconds of day]   [n,2]
%  based on the start of GPS time as a reference:
%     Jan 6, 1980 = [0,0] = [44244,0]
%  This conversion accounts for the leap second 
%  difference between GPS and UTC time.
%
%  Modified Julian Date = JD - 2,400,000.5
%
%  Written: Mike Moreau 12/20/2002

function [mjd,sod] = gps2mjd(week,sec);

% constants
SECONDS_IN_DAY=24*3600;

% Inputs are given as a GPS time in GPS format 
% Must convert to the corresponding UTC time in GPS format
% accounting for leap seconds
[utcweek,utcsec] = utc2gps(gps2utc(week,sec),0);

% Compute total elapsed seconds, using GPS epoch as the reference
total_sec = gpst2sec([utcweek,utcsec]);

% Compute total elapsed days, using GPS epoch as the reference
total_days = floor(total_sec/SECONDS_IN_DAY);

% Compute the MJD corresponding to start of GPS time
% Corresponding to Jan 6, 1980, or total_gps_sec = 0
mjd_gps0 = 44244;  

% Compute the MJD in question
mjd = mjd_gps0 + total_days;

% Compute seconds of day
sod = mod(utcsec,SECONDS_IN_DAY);
