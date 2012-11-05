% Orbit determination of the IBEX mission using the Square Root Information Filter, ESTSRIF
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

% S. Hur-Diaz
% Emergent Space Technlogies, Inc.
% 2/23/2010

% Set up initial conditions
% Osculating Elements at Apogee
GM = JATConstant('muEarth','JGM')*1e-9; %m^3/s^2
RE = JATConstant('rEarth')*1e-3; % km
ra = 48*RE;
rp = 1.4*RE;
OE.sma  = (ra + rp)/2;
OE.ecc = (ra - rp)/(ra + rp);
OE.incl = 45*pi/180;
OE.raan = 0;
OE.argp = 0;
OE.tran = pi;
OE.nu = OE.tran;
xi = kep2cart(OE,GM);

n = sqrt(GM/OE.sma^3);
Tp = 2*pi/n;

% Find perigee time
dtp = mod((2*pi - kepanom(OE,'M')),2*pi)/n;
dtpday = dtp/86400;
oe1 = OE;oe1.tran=1.5*pi;ma1=kepanom(oe1,'M');
dtpday1 = (mod((ma1 - kepanom(OE,'M')),2*pi)/n)/86400;
del_tp = dtpday-dtpday1;
oe1 = OE;oe1.tran=110*pi/180;ma1=kepanom(oe1,'M');
dtpday1 = (mod((ma1 - kepanom(OE,'M')),2*pi)/n)/86400;
del_tp110 = dtpday1-dtpday;

t0 = datenum(2009,5,19); % days

tspan = 0:10*3600:86400*20;

nrev = tspan(end)/Tp
       
% Plot average orbit
[t,x] = integ(@r2bp,tspan,xi,[],GM);

% Initialize estimator options structure
estimOpts = odtbxOptions('estimator');

% Modify default estimator options:
estimOpts = setOdtbxOptions(estimOpts, 'MonteCarloCases', 1);
estimOpts = setOdtbxOptions(estimOpts, 'MonteCarloSeed', 0);
estimOpts = setOdtbxOptions(estimOpts, 'UpdateIterations', 1);
estimOpts = setOdtbxOptions(estimOpts, 'EditRatio', [9 9 9]); 
estimOpts = setOdtbxOptions(estimOpts, 'EditFlag', [2 2 2]); 
estimOpts = setOdtbxOptions(estimOpts, 'refint', 0);

% Initialize measurement options structure
measOptions = odtbxOptions('measurement');

% Modify default measurement options:
gsList  = createGroundStationList; %generates a list of all of the ground stations
gsID = {'AUWS','USAS','USHS'};
gsECEF = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',t0);
end
measOptions = setOdtbxOptions(measOptions,'gsECEF',gsECEF);
measOptions = setOdtbxOptions(measOptions,'epoch',t0); %datenum format
measOptions = setOdtbxOptions(measOptions,'useRange',false); %compute range measurments
measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
measOptions = setOdtbxOptions(measOptions,'useLightTime',false);
measOptions = setOdtbxOptions(measOptions,'useRangerate',true); %compute range rate measurments
measOptions = setOdtbxOptions(measOptions,'useDoppler',false); %don't compute Doppler measurments
measOptions = setOdtbxOptions(measOptions,'rSigma',[.045 .045 .045]*1e-3); % km
measOptions = setOdtbxOptions(measOptions,'useChargedParticle',false);
measOptions = setOdtbxOptions(measOptions,'gsElevationConstraint',0); %elevation constraint
measOptions = setOdtbxOptions(measOptions,'EarthAtmMaskRadius',6478); %km

% Find times of valid obs
[y] = gsmeas(t,x,measOptions); % NB: gsmeas expects states in km
k = any(~isnan(y),1);
tspano = tspan(k);
if tspano(1) ~= tspan(1)
    tspano = [tspan(1) tspano];
end


%% Estimator

Xo = xi;
Po = diag([1e3,1e3,1e3,.1e0,.1e0,.1e0].^2)*1e-6; %km, km/sec
I = eye(length(Xo));
S = I;
C = zeros(0,0,1);

Pnot.Po = Po;
Pnot.Pbaro = S*Po*S';
Xnot.Xo = Xo;
Xnot.Xbaro = S*Xo;
dynarg.tru = GM;
dynarg.est = dynarg.tru;
datarg.tru = measOptions;
datarg.est = measOptions;

Tspan = tspano;

No_pts = length(Tspan)

tic

dynfun.tru = @r2bp;
datfun.tru = @gsmeas;
dynfun.est = dynfun.tru;
datfun.est = datfun.tru;
[t,xhat,Phat,e,y,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,eflg,Pdy,Pdyt]=estsrif(...
    dynfun,datfun,Tspan,Xnot,Pnot,estimOpts,dynarg,datarg,S,C);

elapsed_time = toc;

disp(['Elapsed time is ',num2str(elapsed_time/60),' minutes.'])

% Convert time to days
for i = 1:length(t),
    t{i} = t{i}/86400;
end

plot_results_ric(t,xhat,Phat,e,y,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt,S,1);

nofig = 3;
for i=1:nofig
    figure(i),
    for j=1:3
        subplot(3,1,j),aa=axis;
        axis([t{1}(1) t{1}(end) aa(3) aa(4)])
        hold on,
        tp = dtpday;
        while tp <= Tspan(end)/86400
            plot(tp*[1 1],aa(3:4),'r'),
            plot((tp-del_tp)*[1 1],aa(3:4),'r--')
            plot((tp+del_tp)*[1 1],aa(3:4),'r--')
            tp = tp + Tp/86400;
        end
        hold off
    end
end



