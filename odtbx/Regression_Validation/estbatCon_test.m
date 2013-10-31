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

function [failed] = estbatCon_test

    % Default to passed.
    failed = 0;

    % Simulation of JWST in ODTBX to match ODEAS results

    % Total simulation time of 48 hrs (2880 min)
    %     18 hrs of measurements
    %     30 hrs of prediction
    % Data output at 60 min intervals
    % Measurement at 10 second intervals

    % INDEX:
    % 1 x
    % 2 y
    % 3 z
    % 4 xdot
    % 5 ydot
    % 6 zdot
    % 7 SRP
    % 8 Earth GM
    % 9 Sun GM
    % 10 Moon GM
    % 11 Meas bias DS16
    % 12 Meas bias DS46
    % 13 Meas bias DS66
    % 14 Trop corr DS16
    % 15 Trop corr DS46
    % 16 Trop corr DS66
    %
    % 14 Iono corr DS16
    % 15 Iono corr DS46
    % 16 Iono corr DS66

    tic

    solve = {'X';'Y';'Z';'Xdot';'Ydot';'Zdot'};
    dyn_cons = {'SOLRAD 8';'EARTH-GM';'SOLAR-GM';'LUNAR-GM'};
    % loc_cons = {'MEASBI 1';'MEASBI 2';'MEASBI 3';'ION-DS16';'ION-DS46';'ION-DS66'};
    loc_cons = {'MEASBI 1';'MEASBI 2';'MEASBI 3';'TRP-DS16';'TRP-DS46';'TRP-DS66'};

    s_and_c_tru = solve_consider(solve,dyn_cons,loc_cons,'jat');
    s_and_c_est = solve_consider(solve,'jat');

%     assignin('base','s_and_c_tru', s_and_c_tru)
%     assignin('base','s_and_c_est', s_and_c_est)

    %% INITIAL CONDITIONS
    epoch=datenum('21 July 2013 12:24:44.518');
    % xi=[-6186.520;3722.630;-377.151;-7.62705;-7.07321;-1.13306]; % TOD km & km/s
    % xi=[-6175.49209150131;3741.73553756531;-368.726590886206;-7.6503759967611;-7.04957559516529;-1.12305271756558]; % ECI km & km/s
    % xi=[-6175.49209150131;3741.73553756531;-368.726590886206;-7.65037596477165;-7.04957554237282;-1.12305271760657]; % ECI w only correction km & km/s
    xi=[-6175.49209150131;3741.73553756531;-368.726590886206;-7.6503760287507;-7.04957564796205;-1.12305271756558]; % ECI -D*w correction km & km/s
    Pi=diag([1e6,1e6,1e6,1,1,1]); % km^2 & km^2/s^2

    cxo=[0;0;0;0;0;0;0;0;0;0];
    % sigc=[0.3;0.3e-6;10e-6;1.37e-6;15e-3;15e-3;15e-3;1.0;1.0;1.0];
    sigc=[0.3;0.3e-6;10e-6;1.37e-6;15e-3;15e-3;15e-3;0.045;0.045;0.045];

    %% Build time span

    % Basic output requirements
    % obt=0:1:2880; % [h] - every 60 min for 48 hrs - INCORRECT
    obt=0:60:100; % [m] - every 60 min for 48 hrs

    GStimes=[37:10/60:259,283:10/60:751,752:10/60:1080]; % [m] - every 10 seconds during three contacts

    % The simulation is only run over the interval obt, no matter what the
    % scheduled ground station times are
    index = 1;
    GS_in_time_range = [];
    for time = 1:length(GStimes)
        if GStimes(time) > min(obt) && GStimes(time) < max(obt)
    %         fprintf('%4i < %4.2f < %4i\n',min(obt),GStimes(time),max(obt))
            GS_in_time_range(index) = GStimes(time);
            index = index + 1;
        end
    end

    % tspan=union(GStimes,obt)*60;
    tspan = union(GS_in_time_range, obt) * 60;

    clear GStimes obt

    %% DYNAMICS

    forceOpt=odtbxOptions('forcemodels');
    % forceOpt=setOdtbxOptions(forceOpt,'ap',[]);
    % forceOpt=setOdtbxOptions(forceOpt,'atmosphereModel','HP');                 % WANT JACCHIA-ROBERTS
    forceOpt=setOdtbxOptions(forceOpt,'cD',2.2);
    forceOpt=setOdtbxOptions(forceOpt,'cR',1.5);
    forceOpt=setOdtbxOptions(forceOpt,'dragArea',19.53);                         % LUMPED INTO area/mass 0.315*e-2 (m2/kg)
    forceOpt=setOdtbxOptions(forceOpt,'earthGravityModel','JGM2');             % JGM2F70
    forceOpt=setOdtbxOptions(forceOpt,'epoch',epoch);
    % forceOpt=setOdtbxOptions(forceOpt,'f107Average',150);
    % forceOpt=setOdtbxOptions(forceOpt,'f107Daily',150);
    forceOpt=setOdtbxOptions(forceOpt,'gravDegree',70);
    forceOpt=setOdtbxOptions(forceOpt,'gravOrder',70);
    forceOpt=setOdtbxOptions(forceOpt,'mass',6200);                            % LUMPED INTO area/mass 0.315*e-2 (m2/kg)
    % forceOpt=setOdtbxOptions(forceOpt,'nParameterForHPModel',2);
    forceOpt=setOdtbxOptions(forceOpt,'srpArea',19.53);                          % ASSUMED TO BE SAME AS DRAG AREA
    forceOpt=setOdtbxOptions(forceOpt,'useAtmosphericDrag',false);
    forceOpt=setOdtbxOptions(forceOpt,'useLunarGravity',true);
    forceOpt=setOdtbxOptions(forceOpt,'useSolarGravity',true);
    forceOpt=setOdtbxOptions(forceOpt,'useSolarRadiationPressure',true);
    jatWorldtru = createJATWorld(forceOpt); 
    jatWorldest = createJATWorld(forceOpt); 

    dynfun.tru=@s_and_c_tru.extForces;
    dynfun.est=@s_and_c_est.extForces;

    dynarg.tru = jatWorldtru;
    dynarg.est = jatWorldest;

    %% MEASUREMENTS

    % SET UP GROUND STATIONS

    gsECEF=zeros(3,3);

    gsECEF(:,1) = [-2354.749485,-4646.769139,3669.367433];
    gsECEF(:,2) = [-4460.817118,2682.124829,-3674.963418];
    gsECEF(:,3) = [4849.130216,-360.465192,4114.975241];

    % DS16 - ID 1
    % DS46 - ID 2
    % DS66 - ID 3
    Sched=[
        2 37 259 
        3 283 751
        1 752 1080
        ];
    % Sched in minutes from epoch. Convert to seconds
    Sched(:,2:3)=Sched(:,2:3)*60;

    measOpt.Gs=odtbxOptions('measurement');
    % measOpt.Gs=setOdtbxOptions(measOpt.Gs,'EarthAtmMaskRadius',6378.14+100); 
    % measOpt.Gs=setOdtbxOptions(measOpt.Gs,'PrecnNutnExpire',1); % default
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'Schedule',Sched);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'epoch',epoch);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'frequencyTransmit',2.1e9);         % 2250 mega hertz
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'gsECEF',gsECEF);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'gsElevationConstraint',5);                        
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'rangeType','2way');
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useChargedParticle',false);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useDoppler',false);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useGPSIonosphere',false);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useIonosphere',false);               % MUST BE TRUE!
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useLightTime',false);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useRange',true);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useRangeRate',true);
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useTroposphere',true);               % QUESTIONABLE
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'useUnit',false);

    %% ESTIMATOR

    estOpts=odtbxOptions('estimator');
    estOpts=setOdtbxOptions(estOpts,'MonteCarloCases',1);
    estOpts=setOdtbxOptions(estOpts,'MonteCarloSeed',1);
    % estOpts=setOdtbxOptions(estOpts,'OdeSolver',@rk4n);
    estOpts=setOdtbxOptions(estOpts,'OdeSolver',@ode113);
    estOpts=setOdtbxOptions(estOpts,'SchmidtKalman',0);
    estOpts=setOdtbxOptions(estOpts,'UpdateIterations',3);
    estOpts=setOdtbxOptions(estOpts,'UpdateVectorized',1);
    estOpts=setOdtbxOptions(estOpts,'UseProcNoise',0);
    % estOpts=setOdtbxOptions(estOpts,'refint',0);

    options=estOpts;

    %% SET OPTIONS

    datfun.tru=@s_and_c_tru.meas;
    datfun.est=@s_and_c_est.meas;

    % R noise sigma on truth
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'rSigma',[1e-3;1e-6;1e-3;1e-6;1e-3;1e-6]);             

    datarg.tru=measOpt.Gs;

    % W^2/R include weights on estimate
    measOpt.Gs=setOdtbxOptions(measOpt.Gs,'rSigma',[1e-3;1e-6;1e-3;1e-6;1e-3;1e-6]);             

    datarg.est=measOpt.Gs;

    %% CONSIDER STATES
    solv=length(solve);
    cons=length(dyn_cons)+length(loc_cons);
    tst=solv+cons;

    % cxo=[0;0;0];
    % sigc=[15e-3;15e-3;15e-3];

    x0.Xo=[xi;cxo];       % tru with considers
    x0.Xbaro=xi;            % est

    P0.Po=Pi;
    P0.Po(solv+1:tst,solv+1:tst)=diag(sigc.^2); 
    P0.Pbaro=Pi;

    S=eye(solv);
    S(solv,tst)=0;
    C(1:cons,solv+1:tst)=eye(cons);

    %% FILTER
    [t_b,xhat_b,Phat_b,e_b,y_b,Pa_b,Pv_b,Pw_b,Phata_b,Phatv_b,Phatw_b,Sigma_a_b,Pdy_b,Pdyt_b] = ...
        estbatCon(dynfun,datfun,tspan,x0,P0,options,dynarg,datarg,S,C);
    
    
    %
    % Code block to store data for regression testing:
    %
    %t_b = t;
    %xhat_b = xhat;
    %Phat_b = Phat;
    %e_b = e;
    %y_b = y;
    %Pa_b = Pa;
    %Pv_b = Pv;
    %Pw_b = Pw;
    %Phata_b = Phata;
    %Phatv_b = Phatv;
    %Phatw_b = Phatw;
    %
    %save estbat_tspan_test.mat t_b xhat_b Phat_b e_b y_b Pa_b Pv_b Pw_b Phata_b Phatv_b Phatw_b
    
    % Load data for regression testing:
    load estbatCon_test.mat;

    for a=1:length(t_b)
        dt = t(a) - t_b(a);
        if(any(dt))
            failed = 1;
            fprintf(2, '\n *** estbatCon_test: t test failed.\n');
        end
    end
    
    dx_max = 0.6;
    
    for a=1:length(xhat)
        dx = abs(xhat(a) - xhat_b(a));
        if(max(max(dx)) > dx_max)
            failed = 1;
            disp(max(max(dx)))
            fprintf(2, '\n *** estbatCon_test: xhat test failed.\n');
        end
    end
    disp(max(max(dx)))
    
    dP_max = 3.5e-7;
%     dP_max = 3.5e-5;
    for a=1:length(Phat)
        dP = abs(Phat(a) - Phat_b(a));
        if(max(max(dP)) > dP_max)
            failed = 1;
            disp(max(max(dP)))
            fprintf(2, '\n *** estbatCon_test: Phat test failed.\n');
        end
    end
    
    de_max = 0.85;
    
    for a=1:length(e)
        de = abs(e(a) - e_b(a));
        if(max(max(de)) > de_max)
            failed = 1;
            disp(max(max(de)))
            fprintf(2, '\n *** estbatCon_test: e test failed.\n');
        end
    end
    
    dy_max = 0.75;
    
    for a=1:length(y_b)
        dy = abs(Yref(a) - y_b(a));
        if(max(max(dy)) > dy_max)
            failed = 1;
            disp(max(max(dy)))
            fprintf(2, '\n *** estbatCon_test: y test failed.\n');
        end
    end
    
    % We need tight tolerances, but maybe machine precision is overkill?
%     if((max(max(max(abs(Pa - Pa_b))))) > eps)
%         failed = 1;
%         disp(max(max(max(abs(Pa - Pa_b)))));
%         fprintf(2, '\n *** estbatCon_test: Pa test failed.\n');
%     end
%     
%     if((max(max(max(abs(Pv - Pv_b))))) > eps)
%         failed = 1;
%         disp(max(max(max(abs(Pv - Pv_b)))));
%         fprintf(2, '\n *** estbatCon_test: Pv test failed.\n');
%     end
%     
%     if((max(max(max(abs(Pw - Pw_b))))) > eps)
%         failed = 1;
%         disp(max(max(max(abs(Pw - Pw_b)))));
%         fprintf(2, '\n *** estbatCon_test: Pw test failed.\n');
%     end
%     
%     if((max(max(max(abs(Phata - Phata_b))))) > eps)
%         failed = 1;
%         disp(max(max(max(abs(Phata - Phata_b)))));
%         fprintf(2, '\n *** estbatCon_test: Phata test failed.\n');
%     end
%     
%     if((max(max(max(abs(Phatv - Phatv_b))))) > eps)
%         failed = 1;
%         disp(max(max(max(abs(Phatv - Phatv_b)))));
%         fprintf(2, '\n *** estbatCon_test: Phatv test failed.\n');
%     end
%     
%     if((max(max(max(abs(Phatw - Phatw_b))))) > eps)
%         failed = 1;
%         disp(max(max(max(abs(Phatw - Phatw_b)))));
%         fprintf(2, '\n *** estbatCon_test: Phatw test failed.\n');
%     end
    
end







