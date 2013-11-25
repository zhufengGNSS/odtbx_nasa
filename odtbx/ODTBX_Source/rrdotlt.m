function [y,H,R,t2_lt,x2_lt] = rrdotlt(t1,x1,t2,x2,options)

% RRDOTLT  Calls rrdotang with the light time correction
%
%   [y,H,R] = rrdotlt(t1,x1,t2,x2,options)
% Returns the light time corrected range and range rate (or doppler)
% measurement between two objects and is based on the GEONS DATSIM Math
% Spec.  This function can also include tropospheric and ionospheric delays
% in the correction.
%
% This function can return either '1way' or '2way' measurements based on
% the 'type' specified in the options structure. For '1way' or '1wayFWD'
% measurements, x1 is the spacecraft that receives the range measurements
% at times t. x2 is the transmitting spacecraft/ground station. RRDOTLT
% will interpolate x2's position backward in time to determine the
% measurement. The time of this transmission is calculated with the
% LIGHTTIMECORRECTION function, and reception is t.
%
% For '1wayRTN' measurements, x1 is the spacecraft that sends the range
% measurements at times t. x2 is the receiving spacecraft/ground station.
% RRDOTLT will interpolate x2's position forward in time to determine the
% end point of the measurement. The time of this reception is calculated
% with the LIGHTTIMECORRECTION function, and the time of transmission is t.
%
% For '2way' measurements, x1 is the spacecraft that receives and
% retransmits the range signal from the tracking spacecraft/groundstation.
% x2 is the tracking spacecraft/groundstation that inititates the
% transmission and receives the final measurement. t are the times that x1
% receives the initial transmission. RRDOTLT will interpolate x2's position
% backward in time to determine the time (with the LIGHTTIMECORRECTION
% function) of the initial transmission. RRDOTLT will then interpolate x2's
% position forward in time (with the LIGHTTIMECORRECTION function)to
% determine when it will receive the final measurement.
%
% options is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES             NOTES
%   rangeType           {'1way'(default),'2way',
%                        '1wayFWD','1wayRTN'}    1way=1wayFWD
%   useRange            {true(default), false}
%   useRangeRate        {true(default), false}   can't combine with Doppler
%   useDoppler          {true, false(default)}   can't combine w/ RangeRate
%   frequencyTransmit   {scalar>0, 1.57542e9}    Hz, Only for Doppler and
%                                                measurement errors
%   useGPSIono          {true, false(default)}   only for GPS sats as x2
%   useIono             {true, false(default)}   only for groundstats as x2
%   useTropo            {true, false(default)}   only for groundstats as x2
%   useChargedParticle  {true, false(default)}   only for groundstats as x2
%   epoch                datenum                 UTC time associated with   
%                                                start of simulation.
%
% The measurements are output in y. Each column corresponds to a different
% time. All the measurements for each type are grouped together in rows.
% Thus, using the default options settings, the output y will look like:
%   y = [  range(t1)      range(t2)      range(t3)      range(t4)     ... ;
%          rangeRate(t1)  rangeRate(t2)  rangeRate(t3)  rangeRate(t4) ... ]
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      t1       (1xN)	Times corresponding to x1 (secs)
%      x1       (6xN)   States (km & km/s) of user satellite at times t1
%      t2       (1xN)	Times corresponding to x2 (secs)
%      x2       (6xN)   States (km & km/s) of tracking satellite/groundstat
%                       at times t2
%      options  (1x1)   data structure
%
% OUTPUTS
%      y        (MxN)   measurements
%      H        (Mx6xN) measurement partials matrix
%      R        (MxMxN) measurement covariance
%      t2_lt    (cell)  time outputs from lightTimeCorrection
%      x2_lt    (cell)  state outputs from lightTimeCorrection
%                       Cell dimensions are 1x1 for 1way and 2x1 for 2way
%
% VALIDATION/REGRESSION TEST
%
%  These have been moved to rrdotlt_test.m in the regression testing
%  framework to conform with the new testing format.
%
% keyword: measurement
% See also LOSRANGE, LOSRANGERATE, LOSDOPPLER, RRDOT, GSMEAS, 
%  LIGHTTIMECORRECTION, ODTBXOPTIONS
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
%   Kevin Berry         05/20/2008      Removed light time corrections from
%                                       the rrdot function to help
%                                       simplify, and made rrdotlt to get
%                                       the functionality back into OD
%                                       Toolbox. 
%   Kevin Berry         09/08/2008      Added Validation and Regression
%                                       Tests
%   Kevin Berry         06/25/2009      Added time scale comments
%   Ravi Mathur         08/28/2012      Extracted regression test

type      = getOdtbxOptions(options,'rangeType','1way');
useRange     = getOdtbxOptions(options, 'useRange', true );
useRangeRate = getOdtbxOptions(options, 'useRangeRate', false ); %mutually exclusive with useDoppler
useDoppler   = getOdtbxOptions(options, 'useDoppler', true ); %mutually explusive with useRangeRate
numtypes     = useRange + useRangeRate+useDoppler;

if any( strcmpi( type, {'2way','1wayFWD','1way'} ) )
    direction      = -1; %x2 sends the signal before x1 receives it
    [t2_lt{1}, x2_lt{1}] = lightTimeCorrection(t1, x1, t2, x2, options, direction);
    [y,H,R]        = rrdotang(t1,x1,x2_lt{1},options);  
elseif strcmpi( type, '1wayRTN' )
    direction      = 1; %x2 receives the signal after x1 sends it
    [t2_lt{1}, x2_lt{1}] = lightTimeCorrection(t1, x1, t2, x2, options, direction);
    [y,H,R]        = rrdotang(t1,x1,x2_lt{1},options); 
end
    
if strcmpi( type, '2way' ) 
    direction      = 1; %x2 receives the signal after x1 sends it
    [t2_lt{2}, x2_lt{2}] = lightTimeCorrection(t1, x1, t2, x2, options, direction);
    [y2,H2]            = rrdotang(t1,x1,x2_lt{2},options); 
    % average the pseudorange and rangerate measurements, 
    % Sum Doppler y and H, average others
    if numtypes == 2 && useDoppler
        y = [(y(1,:) + y2(1,:))/2;
               y(2,:) + y2(2,:)      ];
        H = [(H(1,:) + H2(1,:))/2;
               H(2,:)+H2(2,:)        ]; %Sum H for two way Doppler
    elseif numtypes == 1 && useDoppler
        y       = (y+y2);
        H       = H+H2; %Sum H for two way Doppler
    else
        y       = (y+y2)/2;
        H       = (H+H2)/2;
    end
    %rrdotang knows that it is 2way from the options structure, so H was
    %already calculated for a 2way measurement.
end

end