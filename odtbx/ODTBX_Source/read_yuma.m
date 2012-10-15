function [alm] = read_yuma(fn,simdateUTC,toplot)
% READ_YUMA  Read in yuma format GPS almanac file
%
%   [alm] = read_yuma(fn,simdateUTC,toplot)
%
% Read in yuma format GPS almanac file.
% Updates the modulo 1024 GPS week number to the full GPS week
% that corresponds closest to the PC clock.  (I.E. when executed
% in 2002, week 100 would be interpreted as starting 7/21/2001,
% not 12/6/1981.
% Echo summary of almanac data to display.
%
% Input:   fn  - name of yuma formatted almanac file\
%          simdateUTC - the UTC date (in datenum format) of the simulation.
%            This is used when rolling-over the GPS alamanac week number.
%            It will roll-over the week number to make as close to the
%            simulation date as possible.
%            If a date earlier than 1980 is specified no roll-over will be
%            done.
%            If '0' is passed in, this will use the current date as the
%            simulation date (This is backward compatible with previous
%            versions of read_yuma where '0' passed in as a second argument
%            attempted week roll-over to get close to the current date).
%            This argument is optional.  '0' is the default.
%          toplot - 1 - plot almanac data with svmap, 2 plot w/o sv map, 
%            0 - no plots.  This argument is optional.  '0' is the default.
% Output:  alm - [n,13] matrix containing gps almanac data
%          columns of almanac matrix are:
%             PRN
%             Health
%             Eccentricity
%             Time of Applicability  [seconds of week]
%             Inclination [rad]
%             Rate of Right Ascen [rad/s]
%             SQRT(A)  [km^1/2]
%             Longitude of Asc Node at weekly epoch [rad]
%             Argument of Perigee [rad]
%             Mean Anom [rad]
%             Af0 [s]
%             Af1 [s/s]
%             gps week (might be modulo 1023)
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
%   Mike Moreau         05/2002         Original
%   Mike Moreau         06/2003         Updated
%   Kevin Berry         09/2007         Modified for ODTBX
%   Keith Speckman      03/15/2008      Minor comment change
%   Rob Antonucci       03/19/2010      Changed second argument to sim date

if nargin < 1
    fn = input('Enter YUMA GPS almanac file: ','s');
end

% If no simulation date is specified or '0' is specified, add week
% roll-over to given week number based on current date
if (nargin < 2) || (simdateUTC == 0)
    simdateUTC = datenum(clock);
end

if nargin < 3
    toplot = 0;
end

fid = fopen(fn,'r');
i = 0;

line = fgetl(fid);
while feof(fid) == 0
    if length(line)>4  &&  strcmp(line(1:4),'****')
        for count = 1:13
            line = fgetl(fid);
            i = i+1;
            number(i) = str2num(line(27:size(line,2)));
        end
    end
    line = fgetl(fid);
end  % while

fclose(fid);

alm = reshape(number',13,size(number,2)/13)';

% Check for GPS week rollover, a sim date before 1980 means
% don't bother rolling over week
simdateLimit = datenum([1980 01 01 00 00 00]);
if simdateUTC > simdateLimit
    tgps_now = floor(convertTime('GPS','UTC', simdateUTC)/7);
    diff = tgps_now(1) - alm(1,13);
    while diff > 512
        alm(:,13) = alm(:,13) + 1024;
        diff = tgps_now(1) - alm(1,13);
    end
end

% Convert SQRT(A) from m^1/2 to km^1/2
alm(:,7)=alm(:,7)*(1e-3)^(1/2);

if toplot
    disp('yuma plotting function doesn''t work yet')
    disp('there are still a few functions needed that aren''t in ODTBX')
end, return %this line should be removed when it works
toe = alm(:,4);               % time of reference ephemeris
M_not = alm(:,10);            % mean anomaly at toe
ecc = alm(:,3);               % eccentricity
arg_peri = alm(:,9);          % argument of periapsis

%Conversion matrix from PRN to SVN
SVmatrix = [32;61;33;34;35;36;37;38;39;40;46;42;43;41;15;56;...
    0;54;59;51;45;47;60;24;25;26;27;44;29;30;31;0];

% Compute arg of latitude of each PRN
for i = 1:size(alm,1)  % Start loop through each prn

    M_k = M_not(i);

    % Iterate to find eccentric anomaly.
    % Use M_k as initial guess for E.
    % All PRNs iterated the same ammount, some may converge earlier
    % than others, but there is no logic to loof for that.
    E = M_k;                      % First guess for E
    E_old = E+1;                  % Set E-old for at least one loop
    while sum(abs(E - E_old) >= 0.0000000005) > 0
        E_old = E;
        E = M_k + (ecc(i).*sin(E_old));
    end

    % true anomaly as a function of eccentric anomaly
    v = atan2(((sqrt(1 - ecc(i).^2)).*sin(E)),(cos(E) - ecc(i)));

    phi(i) = v + arg_peri(i);		      % argument of latitude
    SVN(i) = SVmatrix(alm(i,1));
end


figure
raan = r2d(gpslan2raan(alm(:,8),alm(:,13)));
ma = r2d(alm(:,10));
ma(find(ma>180)) = ma(find(ma>180))-360;
argp = r2d(alm(:,9));
argp(find(argp>180)) = argp(find(argp>180))-360;
phi = r2d(phi)';
phi(find(phi>180)) = phi(find(phi>180))-360;
phi(find(phi<-180)) = phi(find(phi<-180))+360;

xarg = raan;
%yarg = ma;
%yarg = argp;
yarg = phi;

ind = find(~alm(:,2));
plot(xarg(ind),yarg(ind),'b*')
hold on
ind = find(alm(:,2));
plot(xarg(ind),yarg(ind),'r*')
hold off
ylabel('Arg of Latitude')
xlabel('Right Ascension Node')
title(strrep(sprintf('Almanac: %s',fn),'_','\_'));
axis([0,360,-180,180])
for i = 1:length(SVN)
    if toplot == 1
        t= text(xarg(i)+5,yarg(i),sprintf('%d(%d)',alm(i,1),SVN(i)));
    elseif toplot == 2
        t= text(xarg(i)+5,yarg(i),sprintf('%d',alm(i,1)));
    end
    set(t,'FontSize',8)
end

end
