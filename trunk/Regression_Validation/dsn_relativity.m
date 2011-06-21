function RelLightDelay = dsn_relativity(t3,sc_t3,receive_stn);

% RELLIGHTDELAY   Validation function for relativityLightDelay
%
%  RelLightDelay = dsn_relativity(t3,sc_t3,receive_stn);
%
% Determines upleg and downleg times of a signal sent to a spacecraft.
% Inputs:
% t3: time of reception (Julian days)            		1x1
% sc_t3: spacecraft state in heliocentric (AU & AU/day)       1x6
% receive_stn: station location (ecef m)
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

% Developed by Will Campbell
% Modified (hacked) by Keith Speckman for validation purposes 8/15/08

%function [dt,rec_ID,tran_ID] = lighttime(t3,sc_t3) 
% Receiver Station Determination

% Speed of Light (AU/day)
c = 299792.458/149597870.691*86400;
au = JATConstant('au'); % (m)
mu = [22032; 
    324859;
    398600.44;
    42828;
    126712768;
    37940626;
    5794549;
    6836534;
    982;
    4902.8
    132712440018]/149597870.691^3*86400^2;

    global PM_data
    pm_reader

earth_t3 = jpl(t3,3);
%r23_guess = norm(earth_t3(1:3) - sc_t3(1:3));
%dt23_guess = r23_guess/c;
%sc_t2b_guess = interp1(t_history,SC_history,t3-dt23_guess);
%[site,validity] = stnselect(t3,sc_t2b_guess);
%if validity ~= 1
%    dt = [-1 0];
%    rec_ID = 0;
%    tran_ID = 0;
%    return
%end
%rec_ID = ant(site,I);
%receive_stn = gs_ECEF(rec_ID,:);
%receive_stn = stn_drift(receive_stn,gs_LLA(rec_ID,:),v_site(site,:),stn_epoch,t3);
%% Station Uncertainty
%if drift_flag(5) == 1
%    receive_stn = receive_stn + gs_unc(rec_ID,:).*randoms(1:3);
%end





% ECI Position
stn_ECI = ecef2J2000(t3,receive_stn)/au;
r3_vec = earth_t3(1:3) + stn_ECI;

% 2-3 SOLUTION
% Initialization
delta = 1;
dt23_old = 0;
sc_t2b = sc_t3;

% Iteration
%while delta > 1E-8/86400

    r2_vec = sc_t2b(1:3);
    r23 = norm(r2_vec - r3_vec);
%    ro = lt_delays(t3,sc_t2b,rec_ID,site,M*f1(I));
%    term1 = (r23 + ro)/c;
term1 = 0;
    term2 = 0;
    term3 = 0;
    % Solar Relativistic Effect
%    if relative_flag(11) == 1
%    if 1
        r2 = norm(r2_vec);
        r3 = norm(r3_vec);
        term2 = 2/c^3*mu(11)*log((r2 + r3 + r23 + 2*mu(11)/c/c)/(r2 + r3 - r23 + 2*mu(11)/c/c));
%    end
    % Planetary Relativistic Effects
    for j = 1:10
%        if relative_flag(j) == 1
            rb_t2b = jpl(t3-dt23_old,j);
            rb_t3 = jpl(t3,j);
            r2b = norm(r2_vec - rb_t2b(1:3));
            r3b = norm(r3_vec - rb_t3(1:3));
            term3 = term3 + 2/c/c/c*mu(j)*log((r2b + r3b + r23)/(r2b + r3b - r23));
%        end
    end
    % A glitch was encountered, so real component taken
    dt23 = real(term1 + term2 + term3); 

RelLightDelay = dt23;

end
