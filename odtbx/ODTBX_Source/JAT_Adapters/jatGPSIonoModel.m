function dRho = jatGPSIonoModel(rSc, rGPS, time, tvec )

% JATGPSIONOMODEL  Returns the ionosphere range error for GPS code
% measurement from JAT.
%
%   dRho = jatGPSIonoModel(t_mjd, rSc, rGPS ) returns the range error in
% meters for GPS code measurement. Multiply by -1.0 for carrier phase
% measurement. This model is hard-coded for GPS measurements only and
% expects the time input to be in GPS time.
%
% Acceptable input vector sizes are:
%   rSc(3x1) rGPS(3x1) time(1x1) - 2 satellites at one time
%   rSc(3xN) rGPS(3x1) time(1x1) - multiple sat1 at one time
%   rSc(3x1) rGPS(3xM) time(1x1) - multiple sat2 at one time
%   rSc(3xN) rGPS(3xM) time(1x1) - multiple sat1 and sat2 at one time
%   rSc(3xN) rGPS(3xN) time(1xN) - 2 satellites at multiple times
%
%   INPUTS 
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      rSc          (3xN)	Spacecraft ECI position (m)
%      rGPS         (3xM)   GPS Spacecraft ECI position (m)
%      time         (1xN)   GPS Time in Matlab datenum format (default = 0)
%                           See convertTime for UTC to GPS time conversion
%      tvec         (1x1)   Reference TVEC in electrons/m^2% (default = 2e17)
%
%   OUTPUTS 
%      dRho         (NxM)   ionospheric range error (m)
%
%   keyword: JAT Adapter, measurement model, atmosphere model
%   See also JATTROPOMODEL
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
%   Derek Surka         08/21/2007   	Original
%   Kevin Berry         06/22/2009      Corrected the help to specify the
%                                       input time scale of GPS Time

% Set defaults
if( nargin < 3 )
    time = 0;
end
if( nargin < 4)
    tvec = 2e17;
end

% Check input sizes
n    = size(rSc,2);
m    = size(rGPS,2);
nT   = length(time);

if( nT==1 )
    time = time*ones(1,n);
elseif( nT~=n )
    error('MATLAB:jatGPSIonoModel: The length of the time vector must either equal 1 or the number of initial satellites rSc');
end
tmjd = matlabTime2MJD(time);

dRho = zeros(n,m);

% Get error from JAT
ionoModel = jat.gps.IonoModel(tvec);

for i=1:n 
    iv = ionoModel.Iv(tmjd(i),rSc(:,i));
    for j=1:m 
        dRho(i,j) = ionoModel.error(tmjd(i), rSc(:,i), rGPS(:,j), iv);
    end
end

end
