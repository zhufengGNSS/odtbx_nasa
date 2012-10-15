function [rso_eci,C_1] = getPulsarData(filename,id)
%
% GETPULSARDATA Retrieves pulsar data for use by XNAVMEAS
%
% The following is a sample content of the pulsar data file specified by filename:
%
% % pulsars.txt
% %This is a catalog of pulsars for use with the xnavmeas function
% %
% % Official Name        Nickname                    ECI_Unit (x,y,z)                       Pulse TOA C_1
% %                                       x                  y                  z               (sec^3)
% % =============== ============== ================== ================== ================== =================
%    'B0531+21'          'Crab'     0.104817533369789  0.921144947109562  0.374840327489981 0.000000000207737
%    'B1821-24'          ''         0.094966457808533 -0.902283088049447 -0.420555110432931 0.000050430244604
%    'B1937+21'          ''         0.389831073451460 -0.844226330971047  0.367850018712583 0.000002312225139
%    'J0218+4232'        ''         0.607045484518442  0.417623670844230  0.676081540404973 0.004047715093601
% 
% Any number of comment lines can be had at the beginning of the file
% 
% INPUT: 
%  filename     string variable specifying the pulsar data file name
%  id           a cell array of N strings identifying the pulsar names 
%               Example: id = {'B1937+21','B0531+21'}
% 
% OUTPUT:
%  rso_eci      (3xN) unit vector of each pulsar identified by id
%  C_1          (1xN) time of arrival error parameter
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

% REVISION HISTORY:
% Date         Author            Comment
% 09/23/2010   Sun Hur-Diaz      Original creation
%
%%
fid = fopen(filename);

s = textscan(fid,'%s %s %f %f %f %f','CommentStyle','%');
sl = length(s{1});

fclose(fid);

lid     = length(id);
rso_eci = NaN(3,lid);
C_1     = NaN(1,lid);

for i=1:lid
    % Find the data row with matching name
    ind=strfind(s{1},id{i});
    for j=1:sl
        if ~isempty(ind{j})
            rso_eci(:,i) = [s{3}(j) s{4}(j) s{5}(j)]';
            C_1(1,i)     = s{6}(j);
            break;
        end
    end
    if any(isnan(rso_eci(:,i)))
        error(...
            ['Specified pulsar "',id{i},...
            '" not found in the database "',filename,'"'])
    end
end


    
