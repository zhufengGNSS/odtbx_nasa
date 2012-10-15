function jOptions = jatForcesDemoOptions()

% Force model options for Pre-Release 2 demos
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

% Initialize options structure
jOptions = odtbxOptions('force');

% Choose from any or all of these options and enter desired values
jOptions    = setOdtbxOptions(jOptions, 'epoch', JATConstant('MJDJ2000') );
jOptions    = setOdtbxOptions(jOptions, 'cD', 2.2);
jOptions    = setOdtbxOptions(jOptions, 'cR', 1.2);
jOptions    = setOdtbxOptions(jOptions, 'mass', 1000);
jOptions    = setOdtbxOptions(jOptions, 'dragArea', 20, 'srpArea', 20);
jOptions    = setOdtbxOptions(jOptions, 'earthGravityModel', 'JGM3');
jOptions    = setOdtbxOptions(jOptions, 'gravDegree', 2, 'gravOrder', 2);
jOptions    = setOdtbxOptions(jOptions, 'useSolarGravity', false);
jOptions    = setOdtbxOptions(jOptions, 'useLunarGravity', false);
jOptions    = setOdtbxOptions(jOptions, 'useSolarRadiationPressure', false);
jOptions    = setOdtbxOptions(jOptions, 'useAtmosphericDrag', false);
jOptions    = setOdtbxOptions(jOptions, 'atmosphereModel', 'HP');
jOptions    = setOdtbxOptions(jOptions, 'nParameterForHPModel', 2);
jOptions    = setOdtbxOptions(jOptions, 'f107Daily', 150);
jOptions    = setOdtbxOptions(jOptions, 'f107Average', 150);
jOptions    = setOdtbxOptions(jOptions, 'ap', 15);

return