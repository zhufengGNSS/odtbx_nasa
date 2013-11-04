% GPST2SEC        - convert GPS time to total elapsed seconds since GPS epoch in 1980
%
% [total_gps_seconds] = gpst2sec(GPS_time)
%
% 	Utility function to convert a GPS time structure to pure seconds.  This is
% 	useful when generating time intervals that may cross a week boundary.
%
%	Input:  GPS_time = [gps_weeks gps_secs (from beginning of week)]
%
%	Output:	total_gps_seconds = gps_weeks*86400*7 + gps_secs
%
% 	See also SEC2GPST

function [total_gps_seconds] = gpst2sec(GPS_time)

%	Written by:  Paul M. Stoltz  9/7/97
%  Copyright (c) 1997 by Orion Dynamics and Control, Inc.


%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% verify that there is 1 input variable
if nargin ~= 1     % too many inputs, return with an error message
  fprintf('Only one input is required for GPST2SEC. %d were provided.\n',nargin);
  fprintf('See help on GPST2SEC for details.\n');
  if DEBUG_MODE
    fprintf('Error from GPST2SEC:  ');
    fprintf('Incorrect number of input arguments to GPST2SEC.\n');
    fprintf('Returning to the calling function without any output.\n');
    % return to the calling function without filling in the output variables
    return
  else
    error('Invalid number of inputs to GPST2SEC.\n');
  end % if DEBUG_MODE
end % if nargin ~= 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

total_gps_seconds = GPS_time(:,1).*86400*7 + GPS_time(:,2);

%%%%% END ALGORITHM CODE %%%%%

% end of GPST2SEC


