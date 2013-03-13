function f = rotate_results_test
% POLYDYN_TEST Regression test for polydyn.
%
% F = ROTATE_RESULTS_TEST() runs the regression test
%
% keyword: measurement
% See also rotate_results dcm
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
%   Kenneth Getzandanner    03/13/2013      Original rotate_results_test

load rotate_results_data

[x_out,e_out,P_out]=rotate_results(xhat,e,P,[],'ric');

% Uncomment to generate regression test data
% x_test = x_out;
% e_test = e_out;
% P_test = P_out;
% save ./DataFiles/rotate_results_data x_test e_test P_test xhat e P

% Compare generated and test data
tol = 1e-9;

dx = abs(x_out{1}-x_test{1});
de = abs(e_out{1}-e_test{1});
dP = abs(P_out{1}-P_test{1});

if dx<tol
    f = 0;
else
    f = 1;
    fprintf('ROTATE_RESULTS Regression test failed! dx = %g\n\n', ...
        max(max(dx)))
end

if de>tol
    f = 1;
    fprintf('ROTATE_RESULTS Regression test failed! de = %g\n\n', ...
        max(max(de)))
end

if dP>tol
    f = 1;
    fprintf('ROTATE_RESULTS Regression test failed! dP = %g\n\n', ...
        max(max(dP)))
end

end