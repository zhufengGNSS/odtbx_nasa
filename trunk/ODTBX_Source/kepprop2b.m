function KOEf = kepprop2b(KOE,t,GM)

% KEPPROP2B  Propagates Keplerian elements using 2-body dymanics
%
% INPUTS
%   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
%      KOE       (1x1)	Structure as defined in KEPEL:
%         .sma   (1x1)  Semi-major axis (same length units as GM)
%         .ecc   (1x1)  Eccentricity
%         .incl  (1x1)  Inclination (rad)
%         .raan  (1x1)  Right ascension of the ascending node (rad)
%         .argp  (1x1)  Argument of periapse (rad)
%         .tran  (1x1)  True anomaly (rad)
%      t         (1xN)	Times to propagate KOE to (secs from epoch of KOE)
%      GM        (1x1)  Gravitational constant (same length units as .sma)
% OUTPUTS
%      KOEf      (1x1)  Structure as defined in KEPEL:
%         .sma   (1xN)  Semi-major axis at times t
%         .ecc   (1xN)  Eccentricity at times t
%         .incl  (1xN)  Inclination at times t
%         .raan  (1xN)  Right ascension of the ascending node at times t
%         .argp  (1xN)  Argument of periapse at times t
%         .tran  (1xN)  True anomaly at times t
%
% VALIDATION/REGRESSION TEST
%
%   These tests have been moved to EarthOrbitPlot_test.m to conform to
%   the new regression testing format.
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
%   Kevin Berry         12/15/2008      Original
%   Kevin Berry         06/15/2009      Corrected a bug that let the
%                                       initial mean anomaly be ignored
%   Kevin Berry         06/18/2009      Added self test
%   Ravi Mathur         08/28/2012      Extracted regression test

% These elements are being held fixed over time
KOEf.sma  = KOE.sma*ones(1,length(t));
KOEf.ecc  = KOE.ecc*ones(1,length(t));
KOEf.incl = KOE.incl*ones(1,length(t));
KOEf.raan = KOE.raan*ones(1,length(t));
KOEf.argp = KOE.argp*ones(1,length(t));

% Get the mean anomaly at the epoch of KOE
S0.ecc = KOE.ecc;
S0.nu  = KOE.tran;
M0    = kepanom(S0,'M');

% Propagate true anomaly over times t
T = 2*pi*sqrt(KOE.sma^3/GM); %Orbit period (sec)
M = mod(M0+t/T*2*pi,2*pi); %Mean anomaly (rad) at times t
KOEf.tran = zeros(1,length(t)); %initialize to increase speed

for n=1:length(t)
    S.ecc = KOEf.ecc(n);
    S.M   = M(n);
    KOEf.tran(n) = kepanom(S,'nu');
end

end