function doppler = LOSDoppler(t, r1, v1, r2, v2, options)

% LOSDOPPLER  Line of sight instantaneous doppler between two objects
%
%   doppler = LOSDoppler(r1, v1, r2, v2, f) returns the doppler
% measurement between two vehicles and is simply a conversion of the range 
% rate calculation.
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The parameters that are
% valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   useGPSIono          {true, false(default)} only for GPS sats as x2
%   useIono             {true, false(default)} only for groundstats as x2
%   useTropo            {true, false(default)} only for groundstats as x2
%   useChargedParticle  {true, false(default)} only for groundstats as x2
%   frequencyTransmit   {scalar>0, 1.57542e9}  Hz, Only for Doppler and
%                                              measurement errors
%   epoch                datenum               UTC time associated with 
%                                              start of simulation.
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t         (1xN)	Times corresponding to r1 (secs)
%      r1        (3xN)	User spacecraft position (km)
%      v1        (3xN)  User spacecraft velocity (km/s)
%      r2        (3xN)  Tracking spacecraft/ground station position (km)
%      v2        (3xN)  Tracking spacecraft/ground station velocity (km/s)
%      options   (1x1)  Data structure
%
% OUTPUTS
%      doppler   (1xN)  doppler (Hz)
%
% VALIDATION/REGRESSION TEST
%
%  These have been moved to rrdotlt_test.m in the regression testing
%  framework to conform with the new testing format. 
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, RRDOT, RRDOTLT
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Derek Surka         08/27/2007      Original
%   Kevin Berry         05/08/2008      Replaced with a simpler rdot*f/c
%                                       model
%   Kevin Berry         05/22/2008      Added options structure
%   Kevin Berry         09/08/2008      Added Validation and Regression
%                                       Tests
%   Kevin Berry         09/08/2008      Added Charged Particle Model
%   Kevin Berry         06/25/2009      Added time scale comments
%   Ravi Mathur         08/28/2012      Extracted regression test

c       = JATConstant('c')/1000;
f       = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
rdot    = LOSRangeRate(t, r1, v1, r2, v2, options);
doppler = -rdot*f/c;

end