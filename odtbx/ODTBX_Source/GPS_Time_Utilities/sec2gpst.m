% SEC2GPST        - convert from total seconds since GPS epoch to GPS time
%
% [GPS_time] = sec2gpst(total_gps_secs)
%
% Utility function to convert gps time in seconds from Jan. 6, 1980 to a
% 2-element matrix of gps_week and gps_sec from beginning of week
%
%	Input:  total_gps_seconds = gps_weeks*86400*7 + gps_secs
%			  
%	Output: GPS_time = [gps_weeks gps_secs]
%
%							where gps_secs is from the beginning of the week
%							(week begins at midnight Saturday night)	
%
% See also GPST2SEC

function [GPS_time] = sec2gpst(total_gps_secs)

%	Written by:  Paul M. Stoltz  9/7/97
%  Copyright (c) 1997 by Orion Dynamics and Control, Inc.

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% verify that there is 1 input variable
if nargin ~= 1     % too many inputs, return with an error message
  fprintf('Only one input is required for SEC2GPST. %d were provided.\n',nargin);
  fprintf('See help on SEC2GPST for details.\n');
  if DEBUG_MODE
    fprintf('Error from SEC2GPST:  ');
    fprintf('Incorrect number of input arguments to SEC2GPST.\n');
    fprintf('Returning to the calling function without any output.\n');
    % return to the calling function without filling in the output variables
    return
  else
    error('Invalid number of inputs to SEC2GPST.\n');
  end % if DEBUG_MODE
end % if nargin ~= 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%

GPS_time(:,1) = floor(total_gps_secs./(86400*7));
GPS_time(:,2) = total_gps_secs - GPS_time(:,1).*86400*7;

%%%%% END ALGORITHM CODE %%%%%

% end of SEC2GPST
