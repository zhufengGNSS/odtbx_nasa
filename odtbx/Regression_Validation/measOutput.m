function measOutput(t,y,filename);

% Output measurements from estbat for use in regression testing
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

if( nargin<3 || isempty(filename) )
    filename = 'ODTBX';
end

currentDir=pwd;
stateFileName = strcat(filename,'.meas');
file = fullfile(currentDir,stateFileName);
fid = fopen(file, 'w');

%For Range3D output
% fprintf(fid,'Time (s)\t     Range\n--------    ------------\n');
% x= (y.*2)*10^3;
% for i=1:outputInterval:length(x);
%     if(isnan(x(i)))
%         %Not Printing out to file
%     else
%         fprintf(fid,'%3.4f %17.13f\n', t(i), x(1,i));
%     end
% end
% fclose(fid)

%For rrdot3D output
fprintf(fid,'Time (s)\t     Range (m)\t     RangeRate (m/s)\n--------   --------------   ------------\n');
x(1,:)= (y(1,:).*2)*1e3;
x(2,:)= (y(2,:).*2)*1e3;
for i=1:1:length(x);
    if(isnan(x(1,i)) || x(1,i) == 0)
        %Not Printing out to file
    else
        fprintf(fid,'%3.4f %17.13f %17.13f\n', t(i), x(1,i), x(2,i));
    end
end
fclose(fid);