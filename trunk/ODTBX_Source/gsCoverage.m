function [availTimes PossibleSched] = gsCoverage(dynfun,datfun,x0,t0,dt,tn,options,dynarg,datarg)
% GSCOVERAGE  Checks gsmeas for measurement availability and creates a possible tracking schedule.
%
% AVAILTIMES = GSCOVERAGE(DYNFUN,DATFUN,X0,T0,DT,TN,OPTIONS,DYNARG,DATARG)
% propagates the state X0 over the time T0:DT:TN using DYNFUN, DYNARG, and
% OPTIONS then calls GSMEAS using the propagated state and DATARG. The
% output AVAILTIMES lists all of the time steps where measurements are
% available. The inputs DYNFUN, DATFUN, X0, OPTIONS, DYNARG, and DATARG are
% the same inputs that would be expected for ESTBAT, ESTSEQ, or ESTSRIF.
%
% [AVAILTIMES POSSIBLESCHED] = GSCOVERAGE(...) also provides a possible
% tracking schedule based on all measurement availability over the given
% time interval and step. The schedule is used to restrict the gsmeas
% measurement model to only provide measurements for specific ground
% stations during specific time intervals. The format for each row is:
%      [ gs_index start_time stop_time ].
% The gs_index corresponds to index of the ground stations in gsID or 
% gsECEF. The start_time and stop_time will be in seconds from epoch.
%
% VALIDATION TEST
%
%  To perform a validation test, pass in 'ValidationTest' as the
%  only input argument and specify only one output argument.
%
% REGRESSION TEST
%
%  To perform a regression test, pass in 'RegressionTest' as the
%  only input argument and specify only one output argument.
%
% keyword: measurement
% See also GSMEAS
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
%   Kevin Berry         05/03/2010      Original

%% Determine whether this is an actual call to the program or a test

if strcmpi(dynfun,'ValidationTest')||strcmpi(dynfun,'RegressionTest')
    availTimes = Test_gsCoverage();
else
    [availTimes PossibleSched] = Get_gsCoverage(dynfun,datfun,x0,t0,dt,tn,options,dynarg,datarg);
end
end

%% Main function
function [availTimes PossibleSched] = Get_gsCoverage(dynfun,datfun,x0,t0,dt,tn,options,dynarg,datarg)

%% Check dynfun, datfun, x0, dynarg, and datarg for structured inputs
[dynfun,datfun,x0,dynarg,datarg] = checkInput(dynfun,datfun,x0,dynarg,datarg);

%% Get options from datarg
gsID         = getOdtbxOptions(datarg, 'gsID', [] );
gsECEF       = getOdtbxOptions(datarg, 'gsECEF', []);
useRange     = getOdtbxOptions(datarg, 'useRange', true );
useRangeRate = getOdtbxOptions(datarg, 'useRangeRate', true );
useDoppler   = getOdtbxOptions(datarg, 'useDoppler', false );
numtypes     = useRange + useRangeRate + useDoppler;

datarg       = setOdtbxOptions(datarg,'Schedule',[]); %Clear Schedule

if isempty(gsID)
    for n=1:size(gsECEF,2)
        gsID{n} = num2str(n);
    end
end

%% Create a continuous time span
tspan = t0:dt:tn; 

%% Integrate the state
[tt,xprop] = integ(dynfun,tspan,x0,options,dynarg);

%% Get measurement data
y  = feval(datfun,tspan,xprop,datarg);
y1 = y(1:numtypes:end,:)'*0;
y2 = y1 + repmat(1:length(gsID),size(y1,1),1);

%% Plot the Measurement Times
figure;
plot(tt/3600,y2(:,end:-1:1),'Linestyle','none','Marker','.','MarkerSize',15);
legend(gsID(end:-1:1),'Location','southOutside');
xlabel('Hours');ylabel('Ground Station Number');title('Tracking Coverage');
ylim(gca,[.1 length(gsID)+.9])
set(gca,'YTick',1:length(gsID))

%% Create the array of available times for output
availTimes = tspan(max(~isnan(y))~=0);
PossibleSched = [];
if nargout > 1,
    % Find the times where there is a gap in data at each station
    for n=1:length(gsID)
        numind   = find(~isnan(y1(:,n)));
        nanind   = [0; find(isnan(y1(:,n))); length(y1(:,n))+1];
        strtind  = numind(ismember(numind,nanind+1));
        stopind  = numind(ismember(numind,nanind-1));
        SchedSeg = [n(ones(size(strtind))) tspan(strtind)' tspan(stopind)'];
        
        PossibleSched = [PossibleSched; SchedSeg]; %#ok<AGROW>
    end
end
[~,sortind] = sort(PossibleSched(:,2)); 
PossibleSched = PossibleSched(sortind,:);

end

%% Function to check the inputs for structures
function varargout = checkInput(varargin)
names = {'DYNFUN','DATFUN','X0','DYNARG','DATARG'};
for n=1:length(varargin)
    if isfield(varargin{n},'tru')
        warning('GSCOVERAGE:Struct', [names{n},' was entered as a ',...
            'structure for consider analysis. GSCOVERAGE will continue ',...
            'using only the ".tru" field of ',names{n},'.']);
        varargout{n} = varargin{n}.tru;
    elseif isfield(varargin{n},'est')
        warning('GSCOVERAGE:Struct', [names{n},' was entered as a ',...
            'structure for consider analysis. GSCOVERAGE will continue ',...
            'using only the ".est" field of ',names{n},'.']);
        varargout{n} = varargin{n}.est;
    else
        varargout{n} = varargin{n};
    end
end
end

%% Regression/Validation Test
function failed = Test_gsCoverage()
% Define the initial conditions
x_kep.sma  = 7000; %km
x_kep.ecc  = 0.01; %unitless
x_kep.incl = 35*(pi/180); %radians
x_kep.raan = 250*(pi/180); %radians
x_kep.argp = 270*(pi/180); %radians
x_kep.tran = 100*(pi/180); %radians

muEarth   = 398600.4415; %km^3/s^2
x0 = kep2cart(x_kep,muEarth); %km & km/sec

EpochString = 'April 15, 2010 13:52:24';
epoch = datenum(EpochString);

% Define the dynamics
dynfun = @r2bp;
dynarg = muEarth;

% Define the data
gsList = createGroundStationList();
gsID   = {'HBKS','USHS','USPS','WHSX'};
gsECEF = zeros(3,length(gsID));
for n=1:length(gsID)
    gsECEF(:,n) = getGroundStationInfo(gsList,gsID{n},'ecefPosition',epoch);
end

datfun = @gsmeas;
datarg = odtbxOptions('measurement');
datarg = setOdtbxOptions(datarg,'epoch',epoch);
datarg = setOdtbxOptions(datarg,'gsID',{'HBKS','USHS','USPS','WHSX'});
datarg = setOdtbxOptions(datarg,'rSigma',[1e-2 1e-5 1e-2 1e-5 1e-2 1e-5 1e-2 1e-5]);
datarg = setOdtbxOptions(datarg,'gsECEF',gsECEF);

% Run gsCoverage
[availTimes PossibleSched] = gsCoverage(dynfun,datfun,x0,0,10,60*60*2,[],dynarg,datarg);

% Define expected results
exp_availTimes = [0:10:310 590:10:1180 3080:10:3530 4350:10:4840 5980:10:6560 6770:10:7200];
exp_PossibleSched = ...
          [2           0         310;
           4         590        1180;
           1        3080        3530;
           3        4350        4840;
           2        5980        6560;
           4        6770        7200];
       
passed = all(all(PossibleSched==exp_PossibleSched)) & all(exp_availTimes==availTimes);
failed = ~passed;
if failed
    disp(' ')
    disp('Test Failed!')
else
    disp(' ')
    disp('Test Passed.')
end
 
end