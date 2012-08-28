function failed = varpiles_test()
%
% varpiles_test Regression test for varpiles. Enters demo mode if called
% without an output argument, or regression testing mode if called with
% a single output argument.
%
% See also: varpiles.m
%
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
%
%   Ravi Mathur         08/28/2012      Extracted from varpiles.m

if nargout == 0,
    demomode = true;
else
    demomode = false;
end

% Demo/self-test mode: make up some data
t = 0:100;
lent = length(t);
% Case 1: all delta variances positive:
dVa = ones(1,1,lent);
dVv = ones(1,1,lent);
dVw = ones(1,1,lent);
Va = shiftdim(exp(-t),-1);
Vv = ones(1,1,lent) + 1;
Vw = ones(1,1,lent) + 1;
Vhata = Va - dVa;
Vhatv = Vv - dVv;
Vhatw = Vw - dVw;
V = Va + Vv + Vw;
Vhat = Vhata + Vhatv + Vhatw;
figure % Open new figure window
a = varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat);

if demomode,
    disp('Hit any key to continue')
    pause
else
    for k = length(a):-1:1,
        h1(k) = get(a(k));
    end
end

% Case 2: all delta variances negative:
dVa = -ones(1,1,lent);
dVv = -ones(1,1,lent);
dVw = -ones(1,1,lent);
Va = shiftdim(exp(-t),-1);
Vv = ones(1,1,lent) + 1;
Vw = ones(1,1,lent) + 1;
Vhata = Va - dVa;
Vhatv = Vv - dVv;
Vhatw = Vw - dVw;
V = Va + Vv + Vw;
Vhat = Vhata + Vhatv + Vhatw;
a = varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat);

if demomode,
    disp('Hit any key to continue')
    pause
else
    for k = length(a):-1:1,
        h2(k) = get(a(k));
    end
end

% Case 3: mixed-sign delta variances:
dVa = -ones(1,1,lent);
dVv = ones(1,1,lent);
dVw = zeros(1,1,lent);
Va = shiftdim(exp(-t),-1);
Vv = ones(1,1,lent) + 1;
Vw = ones(1,1,lent) + 1;
Vhata = Va - dVa;
Vhatv = Vv - dVv;
Vhatw = Vw - dVw;
V = Va + Vv + Vw;
Vhat = Vhata + Vhatv + Vhatw;
a = varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat);

if demomode,
    disp('Hit any key to continue')
else
    for k = length(a):-1:1,
        h3(k) = get(a(k));
    end
    close % Close current figure
end

if ~demomode,
    load varpiles_data
    try
        for k = length(h1):-1:1,
            fail1(k) = any(h1(k).YData - h1_save(k).YData);
        end
        for k = length(h2):-1:1,
            fail2(k) = any(h2(k).YData - h2_save(k).YData);
        end
        for k = length(h3):-1:1,
            fail3(k) = any(h3(k).YData - h3_save(k).YData);
        end
        fail = any([fail1,fail2,fail3]);
    catch %#ok<CTCH>
        fail = true;
    end
    failed = fail;
end

end