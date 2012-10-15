function options = initialize_nbody(options)
% INITIALIZE_NBODY Initializes and validates inputs for NBODYPM
% 
%   OPTIONS = initialize_nbody(OPTIONS)
% 
%   Takes in a structure OPTIONS and validates each input given by the
%   user as formatted in ODTBXOPTIONS. After which, the inputs are
%   formatted and initailized for the NBODYPM dynamics function. The
%   function has the option to use Spice files with associated GM values.
% 
%   INPUTS
%   VARIABLE       SIZE    		DESCRIPTION (Optional/Default)
%      options     structure    Options structure from ODTBXOPTIONS
% 
% epoch - Simulation start time (UTC) in DATENUM format. [ scalar >=0 {0} ]
%   Set this property to specify the start time of the simulation. It is
%   used to determine planetary and ECI positions for force models.
% 
% CentralBody - Central body. 
%   [ 'Sun' | 'Jupiter Barycenter' | {'Earth'}]
%   Set this property to specify the central body used in the n-body
%   dynamics model. If given a corresponding spice file and GM, primitive 
%   bodies can be used.
% 
% PointMasses - Point masses. 
%   ['Jupiter Barycenter' | 'Earth Barycenter' | {'Sun'} ]
%   Set this property to specify the point masses used in the n-body
%   dynamics model. If given a corresponding spice file and GM, primitive
%   bodies can be used.
% 
% SpiceFiles - Include Spice files. [ 'which('wld18782.15')' | {[]} ]
%   Set this property to specify the spice files to be included in
%   consideration for point masses and/or the central body. Requires a
%   correlating GM value.
% 
% GM - GM for an associated Spice ID [ '2004660' scalar >= 0 | {[]} ]
%   Set this property to specify the GM for the correlating spice ID. The
%   user is required to provide the spice file in SpiceFiles.
% 
% GM_CB - GM value of the central body [ scalar >= 0 | {3.9860e5} ]
%   Set this property to specify the GM value of the correlating central  
%   body as given in CentralBody. 
% 
% GM_PM - GM value(s) of the point mass(es) [ scalar >= 0 | {1.3271e11} ]
%   Set this property to specify the GM value(s) of the correlating point
%   mass(es) as given in PointMasses.
% 
%   OUTPUTS 
%      options                  Reformatted options sctructure from
%                               ODTBXOPTIONS
% 
% 
%    keyword: initialize, nbody, options 
%    See also ODTBXOPTIONS, SETODTBXOPTIONS, GETODTBXOPTIONS, NBODYPM
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

% Benjamin Asher
% Embry-Riddle Aeronautical University
%
%  REVISION HISTORY
%   Author     			Date         	Comment
%               		(MM/DD/YYYY)
%   Benjamin Asher       05/10/2011     Original

%% Furnish the kernels

cspice_furnsh(which('de421.bsp'))     % Planetary Data
cspice_furnsh(which('naif0009.tls'))  % Leap Seconds information
cspice_furnsh(which('pck00009.tpc')); % PCK Planetary constants kernel
try
    cspice_furnsh(which('de-403-masses.tpc'))
catch %#ok<CTCH>
    % FTP de-403-masses.tpc from the NAIF website
    FTPnaif=ftp('naif.jpl.nasa.gov'); % Connect by annonymous ftp
    cd(FTPnaif,'pub/naif/generic_kernels/pck/'); %go to the folder where the data is held
    fn1='de-403-masses.tpc';
    mget(FTPnaif, fn1);               % download to the current matlab directory
    movefile(fn1,fileparts(which('de421.bsp'))); % move file to the folder where de421.bsp lives
    close(FTPnaif);                   % Closes the ftp connection
    
    cspice_furnsh(which('de-403-masses.tpc'))
end

%% Validate Inputs

epoch = getOdtbxOptions(options,'epoch',[]);
PointMasses = getOdtbxOptions(options,'PointMasses','default');
CentralBody = getOdtbxOptions(options,'CentralBody','default');
GM_Input = getOdtbxOptions(options,'GM',[]);
SpiceFiles = getOdtbxOptions(options,'SpiceFiles',[]);

if isempty(epoch)
    err = MException('NBodyValidation:BadInput', ...
        '\nNo epoch given\n');
    throw(err)
end
if strcmpi('default',PointMasses)
    warnstr = 'No point masses given. Recommend using r2bp. ''Sun'' will be used.\n';
    warning('NbodyValidation:BadInput',warnstr)
    PointMasses = {'Sun'};
end
if strcmpi('default',CentralBody)
    warnstr = 'No central body given. ''Earth'' will be used.\n';
    warning('NbodyValidation:BadInput',warnstr)
    CentralBody = 'Earth';
end
if ~isempty(GM_Input)
    GM_Check = 1;
else
    GM_Check = 0;
end
if ~isempty(SpiceFiles)
    SF_Check = 1;
    cspice_furnsh(SpiceFiles)
else
    SF_Check = 0;
end

%% Initialize Validation

refIds = cspice_spkobj(which('de421.bsp'),999);
refNames = cspice_bodc2s(refIds(:)');
for ii=size(refNames,1):-1:1
    refCell{ii} = deblank(refNames(ii,:));
end

if SF_Check
    refSpiceIDs = cspice_spkobj(SpiceFiles,999);
    refSpice = cspice_bodc2s(refSpiceIDs(:)');
    [row1 ~] = size(refSpice);
    for i=row1:-1:1
        refSpiceCell{i} = deblank(refSpice(i,:));
    end
end

%% Validate PointMasses

for n = length(PointMasses):-1:1
    match1 = strcmpi(PointMasses{n}, refCell);
    indices1 = find(match1 == 1,1);
    if isempty(indices1)
        BadPM(n) = 0;
    else
        BadPM(n) = 1;
        GM_PM{n} = cspice_bodvrd(PointMasses{n}, 'GM', 1);
    end
end

%% Validate CentralBody

match2 = strcmpi(CentralBody, refCell);
indices2 = find(match2 == 1,1);
if isempty(indices2)
    BadCB = 0;
else
    BadCB = 1;
    GM_CB = cspice_bodvrd(CentralBody, 'GM', 1);
end

%% Validate SpiceFiles (if given)

BadPMIndices = find(BadPM == 0);
BadCBIndices = find(BadCB == 0,1);

% Check if Bad inputs are in the SpiceFiles and that the GM is given.
if SF_Check == 1 && (~isempty(BadPMIndices) || ~isempty(BadCBIndices))
    % Check PointMasses names against SpiceFiles
    if ~isempty(BadPMIndices)
        for p = length(BadPMIndices):-1:1
            match3 = strcmpi(PointMasses{BadPMIndices(p)}, refSpiceCell);
            indices3 = find(match3 == 1,1);
            if isempty(indices3)
                BadPM(BadPMIndices(p)) = 0;
            elseif ~isempty(indices3) && GM_Check
                BadPM(BadPMIndices(p)) = 1;
                match4 = strcmpi(PointMasses{BadPMIndices(p)},GM_Input{:,1});
                indices4 = find(match4 == 1,1);
                if ~isempty(indices4)
                    GM_PM{BadPMIndices(p)} = GM_Input{indices4(:,1),2};
                else
                    GM_PM{BadPMIndices(p)} = [];
                end
            elseif ~isempty(indices3) && GM_Check == 0
                err = MException('NbodyValidation:BadInput', ...
                    '\nNo GM given.\n');
                throw(err)
            end
        end
        BadPMIndices = find(BadPM == 0);
    end
    % Check CentralBody names against SpiceFiles
    if ~isempty(BadCBIndices)
        match5 = strcmpi(CentralBody, refSpiceCell);
        indices5 = find(match5 == 1,1);
        if isempty(indices5)
            BadCB = 0;
        elseif ~isempty(indices5) && GM_Check
            BadCB = 1;
            match6 = strcmpi(CentralBody,GM_Input{:,1});
            indices6 = find(match6 == 1,1);
            if ~isempty(indices6)
                GM_CB = GM_Input{indices6(:,1),2};
            else
                GM_CB = [];
            end
        elseif ~isempty(indices3) && GM_Check == 0
            err = MException('NbodyValidation:BadInput', ...
                '\nNo GM given.\n');
            throw(err)
        end
        BadCBIndices = find(BadCB == 0,1);
    end
    
    % Check if GM is given for CentralBody, otherwise the name or ID is invalid
    invalidCB = find(BadCB == 0,1);
    if isempty(invalidCB)
        match7 = isempty(GM_CB);
        indices7 = find(match7 == 1,1);
        if ~isempty(indices7)
            % Display invalid input
            err = MException('NbodyValidation:BadInput', ...
                '\nNo GM given for CentralBody: %s\n',...
                CentralBody);
            throw(err)
        end
    else
        % Display invalid input
        err = MException('NbodyValidation:BadInput', ...
            '\nInvalid Name or ID for CentralBody: %s\n',...
            CentralBody);
        throw(err)
    end
    
    % Check if GM is given for PointMasses, otherwise the name or ID is invalid
    invalidPM = find(BadPM == 0);
    if isempty(invalidPM)
        for i = length(GM_PM):-1:1
            match7(i) = isempty(GM_PM{i});
        end
        indices7 = find(match7 == 1);
        if ~isempty(indices7)
            % Display invalid inputs
            err = MException('NbodyValidation:BadInput', ...
                '\nNo GM given for PointMasses: %s\n',...
                PointMasses{indices7});
            throw(err)
        end
    else
        % Display invalid inputs
        err = MException('NbodyValidation:BadInput', ...
            '\nInvalid Name or ID for PointMasses: %s\n',...
            PointMasses{invalidPM});
        throw(err)
    end
end

% Display invalid inputs
if ~isempty(BadCBIndices)
    err = MException('NbodyValidation:BadInput', ...
        '\nInvalid Name or ID for CentralBody: %s\n',...
        CentralBody);
    throw(err)
elseif ~isempty(BadPMIndices)
    err = MException('NbodyValidation:BadInput', ...
        '\nInvalid Name or ID for PointMasses: %s\n',...
        PointMasses{BadPMIndices});
    throw(err)
end

%% Initialize options for nbodypm.m
options = setOdtbxOptions(options,'epoch',epoch);
options = setOdtbxOptions(options,'PointMasses',PointMasses);
options = setOdtbxOptions(options,'CentralBody',CentralBody);
options = setOdtbxOptions(options,'GM_PM',GM_PM);
options = setOdtbxOptions(options,'GM_CB',GM_CB);

end