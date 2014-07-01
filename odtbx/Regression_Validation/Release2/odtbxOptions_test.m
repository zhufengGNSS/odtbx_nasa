function failed = odtbxOptions_test()
% Regression Test Case
% Function(s) odtbxOptions, setOdtbxOptions, getOdtbxOptions,
% validateOdtbxOptions
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

% test 1a: measurement options struct, set each element
try
    % Creation test
    measOptions = odtbxOptions('measurement');
    measOptions = setOdtbxOptions(measOptions,'rangeType','2way');
    measOptions = setOdtbxOptions(measOptions,'useRange',false);
    
    % get results
    t1a1 = getOdtbxOptions(measOptions,'rangeType','NOT2way');
    t1a2 = getOdtbxOptions(measOptions,'useRange',true);
    
    % check results
    if ~strcmp(t1a1,'2way') || t1a2
        failed = 1;
    end
    
    % check validation
    measOptions = validateOdtbxOptions(measOptions);
catch
    failed = 1;
end

% test 1b: measurement options struct, direct set
try
    measOptionsB = setOdtbxOptions('rangeType','2way','useRange',false);
    
    % get results
    t1b1 = getOdtbxOptions(measOptionsB,'rangeType','NOT2way');
    t1b2 = getOdtbxOptions(measOptionsB,'useRange',true);
    
    % check results
    if ~strcmp(t1b1,'2way') || t1b2
        failed = 1;
    end
catch
    failed = 1;
end

% test 1c: measurement options struct copy - heterogeneous
try
    oldoptsC = odtbxOptions('forceModels');
    measOptionsC = setOdtbxOptions(oldoptsC, measOptionsB);
    
    % get results
    t1c1 = getOdtbxOptions(measOptionsC,'rangeType','NOT2way');
    t1c2 = getOdtbxOptions(measOptionsC,'useRange',true);
    
    % check results
    if ~strcmp(t1c1,'2way') || t1c2 || ~isfield(measOptionsC,'dragArea')
        failed = 1;
    end
catch
    failed = 1;
end

% test 1d: measurement options struct copy - homogeneous
try
    oldoptsD = odtbxOptions('measurement');
    measOptionsD = setOdtbxOptions(oldoptsD, measOptions);
    
    % get results
    t1d1 = getOdtbxOptions(measOptionsD,'rangeType','NOT2way');
    t1d2 = getOdtbxOptions(measOptionsD,'useRange',true);
    
    % check results
    if ~strcmp(t1d1,'2way') || t1d2
        failed = 1;
    end
catch
    failed = 1;
end

% test 2: setting forceModel
try
    forceOptions = odtbxOptions('forceModel');
    forceOptions = setOdtbxOptions(forceOptions,'atmosphereModel','NRL');
    forceOptions = setOdtbxOptions(forceOptions,'earthGravityModel','JGM3');
    forceOptions = setOdtbxOptions(forceOptions,'mass',1500);
catch
    failed = 1;
end

% test 3: setting estimator
try
    estOptions = odtbxOptions('estimator');
    estOptions = setOdtbxOptions(estOptions,'OdeSolver',@jatWorldPropagatorRK8);
    estOptions = setOdtbxOptions(estOptions,'UpdateIterations',5);
catch
    failed = 1;
end

% test 4: incorrectly validates a valid structure
try
    validateOdtbxOptions(measOptions);
    validateOdtbxOptions(forceOptions);
    validateOdtbxOptions(estOptions);
catch
    failed = 1;
end

% test 5: passes an invalid structure
caughtError = false;
try
    validateOdtbxOptions([]);
catch
    caughtError = true;
end
if(~caughtError)
    failed = 1;
end

% test 6: input check
try
    useRange    = getOdtbxOptions(measOptions,'useRange',false);
    mass        = getOdtbxOptions(forceOptions,'mass',200);
    updateIter  = getOdtbxOptions(estOptions,'UpdateIterations',10);
catch
    failed = 1;
end

if( abs(mass-1500)>tol )
    failed = 1;
end
if( abs(updateIter-5)>tol )
    failed = 1;
end
if( useRange )
    failed = 1;
end

% test 7: test single-input creation
op=setOdtbxOptions('MonteCarloCases',5);
if ~isstruct(op)
    failed = 1;
end
if ~isfield(op,'MonteCarloCases')
    failed = 1;
end
if getOdtbxOptions(op,'MonteCarloCases',4) ~= 5
    failed = 1;
end
