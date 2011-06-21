function validationOutput(trajTime,x,measTime,y,outputInterval, validationCase,filename)    	

% validationOutput Output the validation measurements for regression testing.
%   validationOutput(trajTime,x,measTime,y,outputInterval) where the
%   trajTime is the time vector for the trajectory output. x is the
%   trajectory, measTime is the measurement time and y is the matrix which
%   contains the measurements. The outputInterval defines the output the
%   user sets for the trajectory and other items coming out of estbat. 
% 
% validationCase   Measurements   Measurement File   Dynamics File   Main Test Function
%       1       Range & RangeRate     rrdot3D.m         r2bp.m       estbat
% 
% keyword: Estimation, Validation
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

% Kathryn Gregory
% NASA Goddard Space Flight Center



if (validationCase == 1) %First Output type. Output trajectory and measurements (Range & RangeRate)
    trajOutput(trajTime,outputInterval,x,filename);
    measOutput(measTime,y,filename);
% elseif (validationCase == 2) % Next test case
end
