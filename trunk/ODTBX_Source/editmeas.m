function [process,eflags,ME] = editmeas(process,dys,Pdys,eflags,eratios)
% EDITMEAS Measurement editing - scalar or vector
%
%    process = editmeas(dys,Pdys,eflags,eratios,process)
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

% Sun Hur-Diaz
% Emergent Space Technologies, Inc.
% July 19, 2010

% Note that process, dys, Pdys, eflags, and eratios only correspond to measurements
% that the user have 'isel'ected before calling this routine.

% Check for NaNs in the measurement that are selected for processing
% This may or may not have been done outside of this routine, but it
% doesn't hurt to check again
process(isnan(dys)) = false;
ival = find(process == true);
dof = length(ival);

% Perform vector editing
if length(eflags) == 1 && eflags == 1
    
    if length(eratios) ~= 1
        error('EDITMEAS: If Editflag is scalar, so must be Editratio.');
    end
    
    try 
        chistat = chi2inv(erf(sqrt(eratios/2)),dof);
    catch ME
        chistat = chi2inv_odtbx(erf(sqrt(eratios/2)),dof);
    end
    
    if dys(process)' / Pdys(process,process) * dys(process) > chistat
        eflags = 0;
        process(process) = false;
    end
    
else
    
    for k = 1:dof
        ind = ival(k);
        if (eflags(ind)==1) * dys(ind)^2 / Pdys(ind,ind) > eratios(ind)
            eflags(ind) = 0;
            process(ind) = false;
        end
    end

end