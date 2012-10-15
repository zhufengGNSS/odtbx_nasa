% Validation test data for jatStaAzEl.m
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
%   Keith Speckman         05/29/2008   	Original
%   Allen Brown            04/09/2009       Updated data for coincident
%                                           satellite and station (case 4)

Case = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 90*pi/180;
ElEx(Case,1) = 45*pi/180;

Case = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 270*pi/180;
ElEx(Case,1) = 45*pi/180;

Case = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 0*pi/180;
ElEx(Case,1) = 90*pi/180;

Case = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 0*pi/180;
ElEx(Case,1) = 0*pi/180;

Case = 5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 0*pi/180;
ElEx(Case,1) = atan(2/1);

Case = 6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 270*pi/180;
ElEx(Case,1) = -45*pi/180;

Case = 7; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 0*pi/180;
ElEx(Case,1) = 45*pi/180;

Case = 8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case,1) = 0*pi/180;
ElEx(Case,1) = 90*pi/180;

Case = 9; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AzEx(Case+0,1) = 90*pi/180;
ElEx(Case+0,1) = 45*pi/180;
AzEx(Case+1,1) = 270*pi/180;
ElEx(Case+1,1) = 0*pi/180;
AzEx(Case+2,1) = 90*pi/180;
ElEx(Case+2,1) = 0*pi/180;
AzEx(Case+3,1) = 270*pi/180;
ElEx(Case+3,1) = 45*pi/180;
AzEx(Case+4,1) = 45*pi/180;
ElEx(Case+4,1) = -atan(1/sqrt(8));
AzEx(Case+5,1) = 0*pi/180;
ElEx(Case+5,1) = atan(1/2);

