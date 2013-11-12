function doy=gps2doy(GPS_week, GPS_sec)
% converts gps epoch to current day of year

% preallocate UTC
UTC_time(1,1:6)=zeros(1,6);

% compute gpsday and gps seconds since start of GPS time 
gpsday = GPS_week * 7 + GPS_sec ./ 86400;
gpssec = GPS_week * 7 * 86400 + GPS_sec;

% get the integer number of days
total_days = floor(gpsday);

% temp is the number of completed years since the last leap year (0-3)
% the calculation begins by computing the number of full days since
% the beginning of the last leap year.  This is accomplished through
% the rem statement.  Since GPS time started at
% 00:00 on 6 January 1980, five days must be added to the total number
% of days to ensure that the calculation begins at the beginning of a
% leap year.  By subtracting one from this result, the extra day in
% the first year is effectively removed, and the calculation can
% simply be computed by determining the number of times 365 divides
% into the number of days since the last leap year.  On the first day
% of a leap year, the result of this calculation is -1
% so the second statement is used to trap this case.

temp = floor((rem((total_days+5),1461)-1) ./ 365);
I_temp=find(temp < 0);
if any(I_temp), 
  temp(I_temp) = zeros(size(temp(I_temp))); 
end % if

% compute the year
UTC_time(:,1) = 1980 + 4 * floor((total_days + 5) ./ 1461) + temp;

% data matrix with the number of days per month for searching 
% for the month and day
% days in full months for leap year
leapdays =   [0 31 60 91 121 152 182 213 244 274 305 335 366];  
% days in full months for standard year
noleapdays = [0 31 59 90 120 151 181 212 243 273 304 334 365];                                                      

% Leap year flag
% determine which input years are leap years
leap_year = ~rem((UTC_time(:,1)-1980),4);     
I_leap = find(leap_year == 1);                % find leap years
I_no_leap = find(leap_year == 0);             % find standard years

% establish the number of days into the current year
% leap year
if any(I_leap)
  day_of_year(I_leap) = rem((total_days(I_leap) + 5),1461) + 1;   
  doy=day_of_year;
end % if any(I_leap)

% standard year
if any(I_no_leap)
  day_of_year(I_no_leap) = ...
      rem(rem((total_days(I_no_leap) + 5),1461) - 366, 365) + 1;
  doy=day_of_year;
end % if any(I_no_leap)
