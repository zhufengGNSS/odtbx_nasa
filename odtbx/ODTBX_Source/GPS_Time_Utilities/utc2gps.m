function [GPS_week, GPS_sec, GPS_day] = utc2gps(UTC_time,leap_sec)

% UTC2GPS         - Converts UTC time to GPS time
%
% [GPS_week, GPS_sec, GPS_day] = utc2gps(UTC_time,leap_sec);
%                          or
% GPS_time = utc2gps(UTC_time,leap_sec); 
%  
% where GPS_time = [GPS_week GPS_sec]
%
% Converts a UTC time matrix to the equivalent time in GPS weeks,
% GPS seconds, and GPS days (GPS time started at 00:00:00 6 JAN 1980)
%
% Input:  
%   UTC_time - matrix of the form [year month day hour minute second] (nx6)
%               with 4-digit (1980) or 2-digit (80) years,
%               vaild years are 1980 - 2079 (2-digit 80-79)
%   leap_sec - leap seconds applied to UTC relative to GPS (optional)
%               can be a 1x1 or an nx1, if not entered the function will
%               use a look-up table to determine the number of leap seconds
% Output: 
%   GPS_week - GPS week (if 0 or 1 output parameters are used,
%               this is filled with [GPS_week GPS_sec].  See the 
%               alternative calling option from above.
%   GPS_sec  - seconds into the week measured from Sat/Sun midnight
%   GPS_day  - days since the beginning of GPS time (optional)
%
% See also GPS2UTC, GPS2LEAP

% Written by: Maria Evans/Jimmy LaMance 10/9/96
% Copyright (c) 1996 by Orion Dynamics and Control, Inc.

% functions called: UTC2LEAP

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% set minimum and maximum values of acceptable input 
UTC_time_max = [2079 12 31 24 60 60];     % maximum for each UTC time field
UTC_time_min = [1980 1 0 0 0 0];          % minimum for each UTC time field
leap_sec_max = 500;                       % maximum leap second value
leap_sec_max = 0;                         % minimum leap second value

% allocate the output GPS time matrices
GPS_week = ones(size(UTC_time,1),1) * inf;
GPS_sec = ones(size(UTC_time,1),1) * inf;
GPS_day = ones(size(UTC_time,1),1) * inf;

% check the dimensions on the input arguments
% UTC_time input size check
if size(UTC_time,2) ~= 6
  fprintf('\nThe size of UTC_time was %d x %d.\n',size(UTC_time))
  
  if DEBUG_MODE
    fprintf('Error message from UTC2GPS: \n');
    fprintf('The size of the UTC_time input to UTC2GPS must be nx6.\n');
    return
  else
    error('The size of the UTC_time input to UTC2GPS must be nx6.')
  end % if DEBUG_MODE

end % if size(UTC_time,2) ~= 6  

% leap_sec input size check
if nargin >= 2
  if size(leap_sec,2) ~= 1
    fprintf('\nThe size of leap_sec was %d x %d.\n',size(leap_sec))

    if DEBUG_MODE
      fprintf('Error message from UTC2GPS:\n');
      fprintf('The size of the leap_sec input to UTC2GPS must be nx1.\n');
      return
    else
      error('The size of the leap_sec input to UTC2GPS must be nx1.')
    end % if DEBUG_MODE
  end % if size(UTC_time,2) ~= 6  
end % if nargin >= 2  

% UTC_time leap_sec input compatibility check (common n dimension)
if nargin >= 2
  if size(leap_sec,1) ~= 1 & size(leap_sec,1) ~= size(UTC_time,1) 
    fprintf('The size of the UTC_time and leap_sec input to UTC2GPS \n');
    fprintf('must be compatible in the n dimension.\n')
    fprintf('The size of UTC_time and leap_sec were %d x %d and %d x %d.\n',...
             size(UTC_time), size(leap_sec))

    if DEBUG_MODE
      fprintf('Error message from UTC2GPS:\n');
      fprintf('The UTC_time and leap_sec variables in UTC2GPS do not match.\n');
      return
    else
      error('The UTC_time and leap_sec variables in UTC2GPS do not match.')
    end % if DEBUG_MODE

  end % if size(leap_sec,1) ~= 1 | size(leap_sec,1) ~= size(UTC_time,1)  
end % if nargin >= 2

% check the validity (sanity) on the input arguments
% year input (check for 4 or 2 digit year) and convert 2 digit years to 4

% check that all of the years are either 2 or 4 digit (not mixed)
I_year_2_digit = find(UTC_time(:,1) < 1900);

if any(I_year_2_digit)
  % verify that it is really a 2 digit year (0 <= year < 100)
  I_invalid_2_year = find(UTC_time(I_year_2_digit,1) > 100 & ...
                          UTC_time(I_year_2_digit,1) < 0 );
  
  % if there are invalid 2 digit years produce a warning or error
  if any(I_invalid_2_year)
    % if all of the two digit years are invalid, produce an error
    % all invalid (error)
    if length(I_invalid_2_year) == length(I_year_2_digit)    
      if DEBUG_MODE
        fprintf('Error msg: All 2 digit years input to UTC2GPS are invalid.\n');
        return
      else
        error('All 2 digit years input to UTC2GPS are invalid.')
      end % if DEBUG_MODE

    % else if only some are invalid produce a warning 
    % and identify the valid ones
    else
    
      fprintf('Warning: Some 2 digit years input to UTC2GPS are invalid.')
      fprintf('         Invalid inputs will return filled with inf.')

      % convert the valid two digit years to four digit format
      I_valid_2 = find(UTC_time(I_year_2_digit,1) < 100 & ...
                       UTC_time(I_year_2_digit,1) >= 0 ); 
      
      % convert the valid one to 4 digit years 
      % start with the years 1980 - 1999
      I = find(UTC_time(I_year_2_digit(I_valid_2),1) >= 80);
      if any(I)
        UTC_time(I_year_2_digit(I_valid_2(I)),1) = ...
             UTC_time(I_year_2_digit(I_valid_2(I)),1) + 1900;
        clear I
      end % if any(I)

      % finish with the years 2000 - 2079
      I = find(UTC_time(I_year_2_digit(I_valid_2),1) < 80);
      if any(I)
        UTC_time(I_year_2_digit(I_valid_2(I)),1) = ...
             UTC_time(I_year_2_digit(I_valid_2(I)),1) + 2000;
        clear I
      end % if any(I)
    
    end % if length(I_invalid_2_year) == length(I_year_2_digit)
    
  else % if any(I_invalid_2_year)
    % convert the valid one to 4 digit years 
    % start with the years 1980 - 1999
    I = find(UTC_time(I_year_2_digit,1) >= 80);
    if any(I)
      UTC_time(I_year_2_digit(I),1) = UTC_time(I_year_2_digit(I),1) + 1900;
      clear I
    end % if any(I)

    % finish with the years 2000 - 2079
    I = find(UTC_time(I_year_2_digit,1) < 80);
    if any(I)
      UTC_time(I_year_2_digit(I),1) = UTC_time(I_year_2_digit(I),1) + 2000;
      clear I
    end % if any(I)

  end % if any(I_invalid_2_year)
  
end % if any(I_year_2_digit)

% clear out local variables that we're finished with
clear I_year_2_digit I_invalid_2_year I_valid_2

% now all of the valid two digit years are converted to 4 digit years
% complete the checking on the input times

% UTC_time input sanity check
I_year_err = find(UTC_time(:,1) > UTC_time_max(1) | ...
                  UTC_time(:,1) < UTC_time_min(1));
I_month_err = find(UTC_time(:,2) > UTC_time_max(2) | ...
                   UTC_time(:,2) < UTC_time_min(2));
I_day_err = find(UTC_time(:,3) > UTC_time_max(3) | ...
                 UTC_time(:,3) < UTC_time_min(3));
I_hour_err = find(UTC_time(:,4) > UTC_time_max(4) | ...
                  UTC_time(:,4) < UTC_time_min(4));
I_min_err = find(UTC_time(:,5) >= UTC_time_max(5) | ...
                 UTC_time(:,5) < UTC_time_min(5));
I_sec_err = find(UTC_time(:,6) >= UTC_time_max(6) | ...
                 UTC_time(:,6) < UTC_time_min(6));

% print warning messages for invalid inputs
if any(I_year_err)
  fprintf('%d invalid years inputs to UTC2GPS.  First invalid input %d.\n',...
           size(UTC_time(I_year_err,:),1),UTC_time(I_year_err,1));
  fprintf('Invalid entries will return inf for the GPS time.\n');
end % if any(I_year_err)
  
if any(I_month_err)
  fprintf('%d invalid month inputs to UTC2GPS.  First invalid input %d.\n',...
           size(UTC_time(I_month_err,:),1),UTC_time(I_month_err,1));
  fprintf('Invalid entries will return inf for the GPS time.\n');
end % if any(I_month_err)

if any(I_day_err)
  fprintf('%d invalid day inputs to UTC2GPS.  First invalid input %d.\n',...
           size(UTC_time(I_day_err,:),1),UTC_time(I_day_err,1));
  fprintf('Invalid entries will return inf for the GPS time.\n');
end % if any(I_day_err)

if any(I_hour_err)
  fprintf('%d invalid hour inputs to UTC2GPS.  First invalid input %d.\n',...
           size(UTC_time(I_hour_err,:),1),UTC_time(I_hour_err,1));
  fprintf('Invalid entries will return inf for the GPS time.\n');
end % if any(I_hour_err)

if any(I_min_err)
  fprintf('%d invalid minute inputs to UTC2GPS.  First invalid input %d.\n',...
           size(UTC_time(I_min_err,:),1),UTC_time(I_min_err,1));
  fprintf('Invalid entries will return inf for the GPS time.\n');
end % if any(I_min_err)

if any(I_sec_err)
  fprintf('%d invalid second inputs to UTC2GPS.  First invalid input %d.\n',...
           size(UTC_time(I_sec_err,:),1),UTC_time(I_sec_err,1));
  fprintf('Invalid entries will return inf for the GPS time.\n');
end % if any(I_sec_err)

% if any field is invalid, set the entire record to be invalid
if any(I_year_err) | any(I_month_err) | any(I_day_err) | ...
   any(I_hour_err) | any(I_min_err)   | any(I_sec_err)

   
  I_invalid1 = sort([I_year_err; I_month_err; I_day_err; ...
                     I_hour_err; I_min_err;   I_sec_err]);
  
  d_err = diff(I_invalid1);
  I_diff = find(d_err ~= 0);
  
  I_invalid = [I_invalid1(I_diff); I_invalid1(length(I_invalid1))];
  
  % clear out variables we're finished with
  clear I_invalid1 d_err I_year_err I_month_err I_day_err 
  clear I_hour_err I_min_err I_sec_err  
else
  I_invalid = blanks(0);
end  

% if there are any invalid inputs 
if ~isempty(I_invalid)
  % now that we have the indices to all of the invalid inputs, create
  % an index to all of the valid inputs
  % start by creating a matrix of zeros the size of the input
  z1 = zeros(size(UTC_time,1),1);

  % now fill in zeros for all of the valid inputs with the index number for 
  % the invalid ones
  I_invalid_all = z1;
  I_invalid_all(I_invalid) = I_invalid;

  % now all the valid input indices can be established with an or
  v1 = z1 | I_invalid_all;    % sets all the valid indices to 0
  I_valid = find(v1 == 0);

else
  I_valid = 1:size(UTC_time,1);
end % if ~isempty(I_invalid)  

% save the total input so the output can be constructed accordingly
UTC_time_all = UTC_time;
clear UTC_time 

% work with only the valid inputs
if ~any(I_valid)
  if DEBUG_MODE
    fprintf('Error msg: No valid inputs to UTC2GPS.\n')
    fprintf('Returning to calling function with no valid output.\n') 
    
    if nargout == 0 | nargout == 1
      GPS_week = [GPS_week GPS_sec];
    end % if nargout == 0 | nargout == 1
    
    return  
  
  else
    fprintf('Check inputs to UTC2GPS.\n')
    error('No valid inputs to UTC2GPS.')
  end % if DEBUG_MODE
  
end % if ~any(I_valid)

UTC_time = UTC_time_all(I_valid,:);

% compute the number of leap seconds between UTC and GPS at the given UTC time.  
% if there is no leap second input, use the look-up table 
% in the utc2leap function.
if nargin < 2
  % recompute the leap seconds with only the valid times
  clear leap_sec
  leap_sec = utc2leap(UTC_time);
  if isempty(leap_sec)
    fprintf('No data returned from UTC2LEAP.\n');
    fprintf('Verify that the function UTC2LEAP and \n');
    fprintf('the data file LEAPSECS.DAT are in the Matlab path.\n');
    error('Unable to compute the leap second offset between GPS and UTC time');
  end % if isempty(leap_sec)
end % if
%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
% data matrix with the number of days per month             
% days in full months for leap year
leapdays =   [0 31 60 91 121 152 182 213 244 274 305 335];  
% days in full months for standard year 
noleapdays = [0 31 59 90 120 151 181 212 243 273 304 334];                                                     

% Leap year flag  
% determine which input years are leap years
leap_year = ~rem((UTC_time(:,1)-1980),4);     
I_leap = find(leap_year == 1);                % find leap years
I_no_leap = find(leap_year == 0);             % find standard years

% generate a matrix that has the days per completed month for both 
% leap and standard years
if any(I_leap)
  dayspermonth(I_leap) = leapdays(UTC_time(I_leap,2));
end % if any(I_leap)
         
if any(I_no_leap)
  dayspermonth(I_no_leap) = noleapdays(UTC_time(I_no_leap,2));
end % if any(I_no_leap)

% compute the number of leap days encounted in past years, 
% need to add one to the fix computation to get the year correct                                                         
leapyrs = fix((UTC_time(:,1) - 1980) ./ 4) + eval('~leap_year');

% Compute the number of days in completed years since 1980
gpsday = (UTC_time(:,1) - 1980) .* 365 + leapyrs;
                               
% Add the number of days for each completed month
gpsday = gpsday + dayspermonth';

% add the number of days for each completed day 
% (which is 1 less than the current day)
gpsday = gpsday + (UTC_time(:,3) - 1);

% add the fraction of days for each completed hour (hour/24)
% seconds into the current day (avoids round-off errors)
gpssec = UTC_time(:,4) .* 3600; 

% add the fraction of days for each completed minute (minute/1440)
gpssec = gpssec + UTC_time(:,5) .* 60;

% add the fraction of days each completed second (second/86400)
gpssec = gpssec + UTC_time(:,6);

% Subtract 5 days because the starting date is 00:00 6 January 1980
gpsday = gpsday - 5; 

GPS_day(I_valid) = gpsday;        

% formal return parameter GPS week
GPS_week(I_valid) = fix(GPS_day(I_valid) ./ 7); 
days_into_week = floor(rem(GPS_day(I_valid),7));

% Add leap seconds to the gps time 
GPS_sec(I_valid) = days_into_week' * 86400 + gpssec' + leap_sec';

% check to make sure the leap seconds don't force a week rollover
I_next = find(GPS_sec(I_valid) >= 86400 * 7);
if any(I_next)
  GPS_week(I_valid(I_next)) = GPS_week(I_valid(I_next)) + 1;
  GPS_sec(I_valid(I_next)) = GPS_sec(I_valid(I_next)) - 86400 * 7;
end % if

% check the output arguments, if the out requested is length 0 ot 1, 
% return GPS week and GPS seconds in a single vector
if nargout == 0 | nargout == 1
  GPS_week = [GPS_week GPS_sec];
end % if nargout == 0 | nargout == 1

% Return the GPS week, GPS seconds into the week, and the number of 
% completed days (GPS). formal return parameter GPS days since Jan 6 1980
if nargout == 3
  GPS_day(I_valid) = GPS_week(I_valid) * 7 + ...
                     GPS_sec(I_valid) / 86400;	 
end % if nargout == 3                     

%%%%% END ALGORITHM CODE %%%%%

% end UTC2GPS                   
