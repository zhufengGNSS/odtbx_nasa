function y = createJATWorld(jOptions)

% CREATEJATWORLD  Creates a JATModel java object for propagating earth orbits.
%
% CREATEJATWORLD creates a JATModel java object that stores the relevant
% information for propagating Earth-centric orbits using JAT force models.
%
%   y = createJATWorld(jOptions)
%
% Creates a JATModel java object that stores the relevant information for
% propagating Earth-centric orbits using JAT force models. The input is a
% forceModel options data structure created using the SETODTBXOPTIONS 
% function. The output is a java object that is passed to the JATFORCES 
% derivatives function.
%
% Default values are set for the parameters if they do not exist in the
% options structure. Changing the values in this function will change the
% default values for all future runs.
%
%   INPUTS
%   VARIABLE        DESCRIPTION (Optional/Default)
%   jOptions        (Optional) ODTBXOPTIONS structure that contains option name/value 
%                   name/value pairs for use with JAT Force models. Set to
%                   [] (or omit) to use the default ODTBXOPTIONS structure.
%
%   OUTPUTS 
%      y            Java JATModel object
%
%
%   keyword: JAT Forces 
%   See also JATFORCES, ODTBXOPTIONS, SETODTBXOPTIONS, GETODTBXOPTIONS
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
%               		(MM/DD/YYYY)
%   Derek Surka          06/27/2007     Original
%   Derek Surka          09/04/2007     Modified to use OdtbxOptions
%   Kevin Berry          06/22/2009     Added unit comments to the code
%   Ravi Mathur          06/29/2014     Empty input defaults to bare
%                                       ODTBXOptions struct 
                   
% Use default ODTBXOptions structure if needed
if((nargin == 0) || isempty(jOptions))
    jOptions = odtbxOptions();
end

% Get options and Set Defaults
% m            = JATConstant('mjdj2000');
m            = 0;  % this sets the default epoch to 0 Matlab time
epoch        = getOdtbxOptions(jOptions, 'epoch', m); %UTC
mass         = getOdtbxOptions(jOptions, 'mass', 1000); %kg
dragArea     = getOdtbxOptions(jOptions, 'dragArea', 10); %m^2
srpArea      = getOdtbxOptions(jOptions, 'srpArea', 10); %m^2
cD           = getOdtbxOptions(jOptions, 'cD', 2.2);
cR           = getOdtbxOptions(jOptions, 'cR', 0.7);
earthGravity = getOdtbxOptions(jOptions, 'earthGravityModel', '2Body');
degree       = getOdtbxOptions(jOptions, 'gravDegree', 20);
order        = getOdtbxOptions(jOptions, 'gravOrder', 20);

% Check values
if(order>degree)
    warning('ODTBX:createJATWorld', 'Order cannot be greater than degree. Setting order = degree.')
    order = degree;
end

forceFlag(1) = getOdtbxOptions(jOptions, 'useSolarGravity', false);
forceFlag(2) = getOdtbxOptions(jOptions, 'useLunarGravity', false);
forceFlag(3) = getOdtbxOptions(jOptions, 'useAtmosphericDrag', false);
forceFlag(4) = getOdtbxOptions(jOptions, 'useSolarRadiationPressure', false);

% Setup atmospheric model parameters
atmModel     = getOdtbxOptions(jOptions, 'atmosphereModel', 'HP');
atmArray(1)  = getOdtbxOptions(jOptions, 'nParameterForHPModel', 2);
atmArray(2)  = getOdtbxOptions(jOptions, 'f107Daily', 150);
atmArray(3)  = getOdtbxOptions(jOptions, 'f107Average', 150);
atmArray(4)  = getOdtbxOptions(jOptions, 'ap', 4);

y = jat.matlabInterface.ODToolboxJATModel( cR, cD, dragArea, srpArea, mass, matlabTime2MJD(epoch) );
y.initializeForces( forceFlag, earthGravity, degree, order, atmModel, atmArray );

end
