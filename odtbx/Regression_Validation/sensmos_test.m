function fail = sensmos_test
% SENSMOS_TEST Regression test for sensmos
%
% See also sensmos
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Ravi Mathur          08/27/2012       Created: Extracted from sensmos.m

t = (1:12)./24;
for i = length(t):-1:1
    S(:,:,i) = t(i)*ones(6,8);
end
try
    sensmos(t,S)
    fail = false;
catch %#ok<CTCH>
    fail = true;
end

end

