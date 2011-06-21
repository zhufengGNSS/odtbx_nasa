%MATLABJATCOMPARISON This script compares JAT and MATLAB output files to validate the adaptor. 
%   MATLABJATCOMPARISON allows the user to compare JAT output files from the integrator with MATLAB
%   output files from the adaptor which calls the JAT integrator.
%
%    keyword: Plotting Utilities/Graphics, JAT Adaptor, 
% 
%    See also PLOTPOSITIONADAPTOR, PLOTVELOCITYADAPTOR,
%    TESTEOM, CALCULATETRUEANOMALYSTEP, jatRK8
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
%               		(MM/DD/YYYY)
%   Kathryn Bradley     06/06/2006   	Original
%   Emergent Space Technologies

close all
clear all 
clc

[stateFileName, statePathLocation] = uigetfile('*.txt', 'Select MATLAB State File');
stateFile = strcat(statePathLocation, stateFileName);
[stateFileName2, statePathLocation2] = uigetfile('*.txt', 'Select JAT State File');
stateFile2 = strcat(statePathLocation2, stateFileName2);

fid = fopen(stateFile);
stateMatMATLAB = textscan(fid,'%f%f%f%f%f%f%f','headerLines',2);
fclose(fid);
fid2 = fopen(stateFile2);
stateMatJAT = textscan(fid2,'%f%f%f%f%f%f%f','headerLines',0);
fclose(fid2);

epochJat = stateMatJAT{1}(1);
epochsec = 0;

time.JAT = (stateMatJAT{1}-epochJat)*86400; %time in epoch seconds
%time.JAT = stateMatJAT{1};
time.MATLAB = stateMatMATLAB{1};

pos.x1=stateMatMATLAB{2};
pos.y1=stateMatMATLAB{3};
pos.z1=stateMatMATLAB{4};

pos.x2=stateMatJAT{2}(1:length(pos.x1));
pos.y2=stateMatJAT{3}(1:length(pos.y1));
pos.z2=stateMatJAT{4}(1:length(pos.z1));

vel.x1=stateMatMATLAB{5};
vel.y1=stateMatMATLAB{6};
vel.z1=stateMatMATLAB{7};

vel.x2=stateMatJAT{5}(1:length(vel.x1));
vel.y2=stateMatJAT{6}(1:length(vel.y1));
vel.z2=stateMatJAT{7}(1:length(vel.z1));

pos.xdiff = pos.x1-pos.x2;
pos.ydiff = pos.y1-pos.y2;
pos.zdiff = pos.z1-pos.z2;
maxXPosDiff = max(pos.xdiff)
maxYPosDiff = max(pos.ydiff)
maxZPosDiff = max(pos.zdiff)
vel.xdiff = vel.x1-vel.x2;
vel.ydiff = vel.y1-vel.y2;
vel.zdiff = vel.z1-vel.z2;
maxXVelDiff = max(vel.xdiff)
maxYVelDiff = max(vel.ydiff)
maxZVelDiff = max(vel.zdiff)
pos.totdiff = sqrt(pos.xdiff.^2+pos.ydiff.^2+pos.zdiff.^2);
vel.totdiff = sqrt(vel.xdiff.^2+vel.ydiff.^2+vel.zdiff.^2);

type = 'State: ';

plotpositionadaptor(time,pos,type);
plotvelocityadaptor(time,vel,type);
