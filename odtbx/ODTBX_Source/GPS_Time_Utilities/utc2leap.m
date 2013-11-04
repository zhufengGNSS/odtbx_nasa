function [leap_sec] = utc2leap(UTC_time)

% UTC2LEAP        - Determines the number of leap seconds of offset between GPS and UTC time
%
% [leap_sec] = utc2leap(UTC_time);
%
% Determines the number of leap seconds of offset between GPS and UTC time.
%
% Input:  
%   UTC_time - matrix of the form [year month day hour minute second]
%               with 4-digit year (1980), nx6 matrix
% Output: 
%   leap_sec - leap seconds relating UTC to GPS time
%
% See also UTC2GPS, GPS2UTC

% Written by: Jimmy LaMance 10/9/96
% Copyright (c) 1996 by Orion Dynamics and Control, Inc.

% functions called: none

% data files called: leapsecs.dat
%   leapsecs.dat has the form of [year month day hour minute second leapsec]
%   where the year month day hour minute second are UTC times, 
%   year is four digits

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug variable
global DEBUG_MODE

% set minimum and maximum values of acceptable input 
UTC_time_max = [2079 12 31 24 60 60];     % maximum for each UTC time field
UTC_time_min = [1980 1 0 0 0 0];          % minimum for each UTC time field

% check the dimensions on the input arguments
% UTC_time input size check
if size(UTC_time,2) ~= 6
  fprintf('\nThe size of UTC_time was %d x %d.\n',size(UTC_time))
  
  if DEBUG_MODE
    fprintf('Error message from UTC2LEAP: \n');
    fprintf('The size of the UTC_time input to UTC2LEAP must be nx6.\n');
    return
  else
    error('The size of the UTC_time input to UTC2LEAP must be nx6.')
  end % if DEBUG_MODE

end % if size(UTC_time,2) ~= 6  

% allocate the leap_sec output
leap_sec = ones(size(UTC_time,1),1) * inf;

% check the validity (sanity) on the input arguments
% year input (check for 4 or 2 digit year) and convert 2 digit years to 4

% check that all of the years are either 2 or 4 digit (not mixed)
I_year_2_digit = find(UTC_time(:,1) < 1900);
if any(I_year_2_digit)
  % verify that it is really a 2 digit year (0 <= year < 100)
  I_invalid_2_year = find(UTC_time(I_year_2_digit,1) > 100 & UTC_time(I_year_2_digit,1) < 0 );
  
  % if there are invalid 2 digit years produce a warning or error
  if any(I_invalid_2_year)
    % if all of the two digit years are invalid prodice an error
    if length(I_invalid_2_year) == length(I_year_2_digit)    % all invalid (error)
      if DEBUG_MODE
        fprintf('Error msg: All 2 digit years input to UTC2LEAP are invalid.\n');
        return
      else
        error('All 2 digit years input to UTC2LEAP are invalid.')
      end % if DEBUG_MODE

    % else if only some are invalid produce a warning 
    % and identify the valid ones
    else
    
      fprintf('Warning: Some 2 digit years input to UTC2LEAP are invalid.')
      fprintf('         Invalid inputs will return filled with inf.')

      % convert the valid two digit years to four digit format
      I_valid_2 = find(UTC_time(I_year_2_digit,1) < 100 & UTC_time(I_year_2_digit,1) >= 0 ); 
      
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
I_year_err = find(UTC_time(:,1) > UTC_time_max(1) | UTC_time(:,1) < UTC_time_min(1));
I_month_err = find(UTC_time(:,2) > UTC_time_max(2) | UTC_time(:,2) < UTC_time_min(2));
I_day_err = find(UTC_time(:,3) > UTC_time_max(3) | UTC_time(:,3) < UTC_time_min(3));
I_hour_err = find(UTC_time(:,4) > UTC_time_max(4) | UTC_time(:,4) < UTC_time_min(4));
I_min_err = find(UTC_time(:,5) >= UTC_time_max(5) | UTC_time(:,5) < UTC_time_min(5));
I_sec_err = find(UTC_time(:,6) >= UTC_time_max(6) | UTC_time(:,6) < UTC_time_min(6));

% if any field is invalid, set the entire record to be invalid
if any(I_year_err) | any(I_month_err) | any(I_day_err) | ...
   any(I_hour_err) | any(I_min_err)   | any(I_sec_err)

   
  I_invalid1 = sort([I_year_err; I_month_err; I_day_err; ...
                     I_hour_err; I_min_err;   I_sec_err]);
  
  d_err = diff(I_invalid1);
  I_diff = find(d_err ~= 0);
  
  I_invalid = [I_invalid1(I_diff); I_invalid1(length(I_invalid1))];
  
  % clear out variables we're finished with
  clear I_invalid1 d_err I_year_err I_month_err I_day_err I_hour_err I_min_err I_sec_err
end  

% if there are any invalid inputs 
if exist('I_invalid')
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
end % if  

% save the total input so the output can be constructed accordingly
UTC_time_all = UTC_time;
clear UTC_time 

% work with only the valid inputs
if ~any(I_valid)
  if DEBUG_MODE
    fprintf('Warning: There are no valid input times to UTC2LEAP.\n')
    fprintf('Returning to calling function with no valid output.\n') 
    return  
  
  else
    fprintf('Warning: There are no valid input times to UTC2LEAP.\n')
    fprintf('Returning to calling function with no valid output.\n') 
    error('No valid inputs to UTC2LEAP.')
  end % if DEBUG_MODE
  
end % if ~any(I_valid)

UTC_time = UTC_time_all(I_valid,:);

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
% verify that the leapsecs.dat file exist
if exist('leapsecs.dat') ~= 2
  fprintf('Unable to open the data file LEAPSECS.DAT.\n');
  fprintf('Verify that the file is in the Matlab path.\n');
  if DEBUG_MODE == 1
    fprintf('Error msg: Unable to open LEAPSECS.DAT.\n')
  else 
    error('Unable to open LEAPSECS.DAT.');
  end % if DEBUG_MODE == 1
end % if exist('leapsecs.dat') ~= 2

% load leap seconds data file (leapsecs.dat)
load leapsecs.dat

% verify that the leapsecs loaded has the correct number of columns
if size(leapsecs,2) ~= 7
  fprintf('Error in LEAPSECS.DAT file.\n');
  fprintf('There were %d columns in the file.  There should be 7 columns.\n',size(leapsecs,2));
  fprintf('Verify the Matlab path is correct and the correct version of LEAPSECS.DAT is in use.\n');
  if DEBUG_MODE == 1
    fprintf('Error msg: Invalid LEAPSECS.DAT file.\n')
  else 
    error('Invalid LEAPSECS.DAT file.');
  end % if DEBUG_MODE == 1
end % if size(leapsecs,2) ~= 7
    
n_leaps = size(leapsecs,1);   % number of leap seconds given in the data file

% search leapsecs for each of the time given in UTC_time. Only have to search based on 
% year and month since the leap seconds are added at the beginning of Jan and July.                     
for jjj = n_leaps:-1:1
  I_n = find(UTC_time(:,1) + UTC_time(:,2) / 12 < leapsecs(jjj,1) + leapsecs(jjj,2) / 12 );
  leap_sec(I_n) = ones(size(I_n,1),1) * leapsecs(jjj,7) - 1;
end % for jjj 

% find 
I_l = find(UTC_time(:,1) + UTC_time(:,2) / 12 >= leapsecs(n_leaps,1) + leapsecs(n_leaps,2) / 12 );
leap_sec(I_l) = ones(size(I_l,1),1) * leapsecs(n_leaps,7);

%%%%% END ALGORITHM CODE %%%%%

% end of UTC2LEAP
