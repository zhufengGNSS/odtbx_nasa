function [t,x] = read_spephem(fn)

%READ_SPEPHEM Read in SP Ephemeris formatted ASCII text file
% 
%  [t,x] = read_spephem(fn)
%
%  Input:
%     fn          name of SP Ephemeris data file
%  Outputs:
%     tSim        Matlab datenum time (1xN)
%     x           ECI state of the satellite (6xN)
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
%   Kevin Berry         01/2008         Original

%% open file for reading
fid = fopen(fn,'r');

%% read in the data
% Format should be [time X Y Z Vx Vy Vz]
rawdata = fscanf(fid, '%f %f %f %f %f %f %f', [7 inf]);

%% Convert the time column into Matlab time
for n=1:size(rawdata,2)
    YYDDDHHMMSS(n,:)=num2str(rawdata(1,n), '%015.3f');
end
t=datenum(YYDDDHHMMSS(:,[1:2 6:end]), 'yyHHMMSS.FFF')'+str2num(YYDDDHHMMSS(:,3:5))'-1;
%% Pull out position and velocity
x=rawdata(2:7,:);

%% Close file
fclose(fid);





