function qt = slerp(q0, q1, t)
% Interpolates between two quaterntions (uniform angular velocity, fixed axis).
%
% Slerp stands for spherical linear interpolation.  It refers to constant
% speed motion along a unit radius great circle arc, where the arc can
% exist in any number of dimensions.  Here it is applied to quaternion
% rotations.
%
% Inputs:
%   q0      4x1        Initial Quaternion of the form: [vector;scalar]
%   q1      4x1        Final Quaternion of the form: [vector;scalar]
%   t       1x1        Linear parameter of the "distance" to interpolate
%                      from q0 to q1, 0<=t<=1, where t=0 is exactly q0 and
%                      t=1 is exactly q1.
% Outputs:
%   qt      4x1        Resulting Quaternion of the form: [vector;scalar]
%
% References:
% http://en.wikipedia.org/wiki/Slerp
% http://number-none.com/product/Understanding%20Slerp,%20Then%20Not%20Using%20It/
% http://www.alecjacobson.com/weblog/?p=981
% http://theory.org/software/qfa/writeup/node12.html
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

% we must rigidly enforce the unit length assumption here
q0 = q0/norm(q0);
q1 = q1/norm(q1);

% the cosine of the angle between the two vectors (in 4D space)
d = dot(q0,q1);

if (d > 0.9995)
    % if the angle is is too close, then inearly interpolate between the
    % two
    qt = q0 + t*(q1 - q0);
else
    % clamp d to -1 to 1 for acos
    d = max(-1,d);
    d = min(d,1);
    
    theta_0 = acos(d); % the angle between the two vectors (in 4D space)

    % the linear interpolation of the angle between the first vector and
    % the result
    theta = theta_0*t;
    
    q2 = q1 - q0*d;
    q2 = q2/norm(q2);
    
    qt = q0*cos(theta) + q2*sin(theta);
end    
qt = qt/norm(qt); % always normalize
