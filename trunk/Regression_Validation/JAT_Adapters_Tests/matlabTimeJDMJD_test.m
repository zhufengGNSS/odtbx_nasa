function failed = matlabTimeJDMJD_test()
% Regression Test Case
% Function(s) matlabTime2JD, matlabTime2MJD, jD2MatlabTime, mJD2matlabTime,
% jD2MJD, mJD2JD
%
% This function tests the various conversions to/from Matlab time, JD, MJD
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

failed  = 0;
tol     = 1e-10;
epoch   = datenum('June 1, 2007');
jD      = 2.454252500000000e+006;
mJD     = 54252;

outJD   = matlabTime2JD(epoch);
outMJD  = matlabTime2MJD(epoch);

if( abs( jD-outJD ) > tol )
    failed = 1;
end

if( abs( mJD-outMJD ) > tol )
    failed = 1;
end

if( abs( epoch-jD2MatlabTime(outJD) ) > tol )
    failed = 1;
end

if( abs( epoch-mJD2MatlabTime(outMJD) ) > tol )
    failed = 1;
end

if( abs( mJD-jD2MJD(outJD) ) > tol )
    failed = 1;
end

if( abs( jD-mJD2JD(outMJD) ) > tol )
    failed = 1;
end


