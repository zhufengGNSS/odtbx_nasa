function y = convertTime(to, from, time)

% CONVERTTIME  Converts the input time from one timescale to another
% as specified by the user.
%
%   y = convertTime(to, from, time) converts the input time from one
% timescale to another as specified by the user. The input and output times
% are in datenum format. Valid input and output timescales can be found by
% calling convertTime with no input arguments.  Valid time scales are:
% 'GPS', 'UTC', 'TT', 'TDB', and 'TAI'.
%
% Special handling for GPS times:
% Note that GPS time inputs and outputs are based on the GPS epoch of
% Jan 6, 1980 and not the standard datenum epoch.  Calling this fuction
% with GPS inputs or outputs will require a correction to be applied to
% any datenum-based times, as in:
%
% for GPS time input that starts in datenum format: 
% nowdatenum = now;
% nowgpsinput = nowdatenum - datenum('Jan 6 1980');
% nowutcoutput = convertTime('UTC','GPS',nowgpsinput);
%
% for GPS time output that needs to be converted to datenum format:
% nowdatenum = now;
% nowgpsoutput = convertTime('GPS','UTC',nowdatenum);
% nowutcdatenum = nowgpsoutput + datenum('Jan 6 1980');
%
% Note that datenum('Jan 6 1980') = 723186 (exact).
%
%   INPUTS
%   VARIABLE    SIZE    	DESCRIPTION (Optional/Default)
%     to        (1X1)       String specifiying output timescale
%     from      (1X1)       String specifying input timescale
%     time      (1XN)       Input time in datenum format
%
%   OUTPUTS
%      y        (1XN)       Output time in datenum format
%
%   keyword: JAT Time
%   See also DATENUM, MATLABTIME2MJD, JD2MATLABTIME, MATLABTIME2JD, JD2MJD,
%   MJD2JD
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2011 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

%   REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Derek Surka         08/07/2007   	Original
%   Derek Surka         09/24/2007      Accept input time vector
%   Kevin Berry         10/2008         Added GPS origin date to epoch (Jan 6, 1980)
%   Keith Speckman      03/18/2008      Corrected GPS epoch modification so that time
%                                         was adjusted BEFORE mjd was calculated

validTimeScales = {  % IMPORTANT: do not switch the order of these entries
    'GPS'
    'UTC'
    'TT'
    'TDB'
    'TAI'};

if( nargin == 0 )
    disp('Valid input and output timescales are:');
    disp(validTimeScales)
    return
end

if( nargin < 3 )
    error('MATLAB:convertTime: There must be 3 input arguments to convertTime.');
end

if( ~isnumeric(time) )
    error('MATLAB:convertTime: Input time must be in datenum format.');
end

% the delta between GPS epoch of Jan 6 1980 and the datenum epoch, faster
% than always calling datenum.
GPSDatenumBias = 723186; % = datenum('Jan 6 1980'); % exact

toIndex     = getIndex(to,validTimeScales);
fromIndex   = getIndex(from,validTimeScales);

if ( fromIndex == 1 )
    time    = time + GPSDatenumBias; % Corrects for GPS Epoch. 
end

mjd         = matlabTime2MJD(time);   % JAT time conversions use MJD

if( toIndex == fromIndex )
    y = time;
else
    yMJD = convert(toIndex,fromIndex,mjd);
    y    = mJD2MatlabTime(yMJD);
end


if ( toIndex == 1 )
    y    = y - GPSDatenumBias; % Corrects for GPS Epoch. 
end

end

%-------------------------------------------------------------------------
function outputTime = convert(toIndex,fromIndex,mjd)

nT          = length(mjd);
outputTime  = zeros(1,nT);

switch( toIndex )
    case 1 % GPS
        switch( fromIndex )
            case 5
                outputTime = mjd - jat.spacetime.TimeUtils.TAI_GPS/86400;
            otherwise
                tai = convert(5,fromIndex,mjd);
                outputTime = convert(toIndex,5,tai);
        end
    case 2 % UTC
        switch( fromIndex )
            case 1
                for k=1:nT
                    outputTime(k) = jat.spacetime.TimeUtils.gps2utc(mjd(k));
                end
            case 3
                for k=1:nT
                    outputTime(k) = jat.spacetime.TimeUtils.TTtoUTC(mjd(k));
                end
            case 4
                tt = convert(3,fromIndex,mjd);
                outputTime = convert(toIndex,3,tt);
            case 5
                for k=1:nT
                    outputTime(k) = mjd(k) - jat.spacetime.TimeUtils.tai_utc(mjd(k))/86400;
                end
        end
    case 3 % TT
        switch( fromIndex )
            case 1
                outputTime = mjd + ( jat.spacetime.TimeUtils.TAI_GPS + jat.spacetime.TimeUtils.TT_TAI )/86400;
            case 2
                for k=1:nT
                    outputTime(k) = jat.spacetime.TimeUtils.UTCtoTT(mjd(k));
                end
            case 4
                for k=1:nT
                    outputTime(k) = jat.spacetime.TimeUtils.TDBtoTT(mjd(k));
                end
            case 5
                outputTime = mjd + jat.spacetime.TimeUtils.TT_TAI/86400;
        end

    case 4 % TDB
        switch( fromIndex )
            case 3
                for k=1:nT
                    outputTime(k) = jat.spacetime.TimeUtils.TTtoTDB(mjd(k));
                end
            otherwise
                tt = convert(3,fromIndex,mjd);
                outputTime = convert(4,3,tt);
        end

    case 5 % TAI
        switch( fromIndex )
            case 1
                outputTime = mjd + jat.spacetime.TimeUtils.TAI_GPS/86400;
            case 2
                for k=1:nT
                    outputTime(k) = mjd(k) + jat.spacetime.TimeUtils.tai_utc(mjd(k))/86400;
                end
            case 3
                outputTime = mjd - jat.spacetime.TimeUtils.TT_TAI/86400;
            case 4
                tt = convert(3,fromIndex,mjd);
                outputTime = convert(5,3,tt);
        end

    otherwise
        error('MATLAB:convertTime: Index out of bounds.')
end

end

