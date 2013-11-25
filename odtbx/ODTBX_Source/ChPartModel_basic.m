function [ drChPart ] = ChPartModel_basic( t, r1, r2, options )
% Calculate the signal delay (in km) caused by a signal traveling 
% through the Charged Particles above the Ionosphere.  This function will 
% only be accurate if r2 is a ground station.
%
% [ drChPart ] = ChPartModel_basic( t, r1, r2, options ) 
%
% INPUTS
%   t         (1xN)	Times corresponding to r1 (secs)
%   r1        (3xN)	User spacecraft position (km) in ECI
%   r2        (3xN) Tracking spacecraft/ground station position (km) ECI
%   options   (1x1) Data structure
% OUTPUTS
%   drChPart  (1xN) Range error in meters due to Charged Particles
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
%   Jason Schmidt       11/7/2013       original ChartModelModel_basic


    height_ionosphere = 3000; %assumption is the Ionosphere extends out to 3000 km altitude
    % If r2 is not below the Ionosphere, set drChPart to 0
    Re = JATConstant('rEarth')/1000; % Radius of the earth in km
    if all((sqrt((sum(r2.^2)))-Re) > height_ionosphere) %
        drChPart = zeros(1,numel(t));
        return
    end
    
    
    % If r1 is not above the Ionosphere, set drChPart to 0
    Re = JATConstant('rEarth')/1000; % Radius of the earth in km
    if all((sqrt((sum(r1.^2)))-Re) < height_ionosphere)
        drChPart = zeros(1,numel(t));
        return
    end
    %r1 and r2 states must be in ECI
    epoch                     = getOdtbxOptions(options, 'epoch', datenum('Jan 1 2006') ); %UTC
    Ephem.SignalFreq          = getOdtbxOptions(options,'frequencyTransmit',JATConstant('L1Frequency'));
    Ephem.Epoch               = epoch+t/86400; %UTC
    Ephem.satPos              = r1*1000;
    drChPart                  = (jatChargedParticleModel(Ephem)/1000)';
    drChPart(isnan(drChPart)) = 0;
    
end

