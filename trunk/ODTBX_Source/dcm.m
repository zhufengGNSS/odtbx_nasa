function D = dcm(type, varargin)
% DCM  Direction Cosine Matrix.
%   D = DCM(TYPE,...) returns a 3x3 direction cosine matrix, or a 3x3xN
%   array of DCMs, of the TYPE the user specifies. The first 3 characters
%   of TYPE, regardless of case, must be unambiguous. Other input
%   arguments depend on the TYPE, as specified below:
%      TYPE  INPUTS                 DESCRIPTION
%      'ax1' rotation angle         Rotation (rad.) about 1st Basis Vector
%      'ax2' rotation angle         Rotation (rad.) about 2nd Basis Vector
%      'ax3' rotation angle         Rotation (rad.) about 3rd Basis Vector
%      'axs' rot. angle and axis    Rot. (rad.) about Specified Axis Vector
%                                   Note that the user must ensure that the
%                                   magnitude of the specified vector is
%                                   unity. Also, the size of the vector
%                                   must be 3x1.
%      'ric' position and velocity  Radial, Intrack, Crosstrack
%                                   Note that the position and velocity
%                                   vectors must both be specified as
%                                   3x1 vectors.
%      'vnb' position and velocity  Velocity, Normal, Bi-normal
%                                   Note that the position and velocity
%                                   vectors must both be specified as
%                                   3x1 vectors.
%      'los' pos. 1 and pos. 2      Line-of-Sight, Normal, Bi-normal
%                                   Note that both position vectors
%                                   must be specified as 3x1 vectors.
%      'con' relative pos. & vel.   Conjunction Basis
%                                   Note that the relative position and
%                                   velocity vectors must both be specified
%                                   as 3x1 vectors.
%   In all cases, the DCM has as its rows the "new" basis vectors
%   corresponding to TYPE, i.e. the DCM is "from" the original basis to the
%   "new" TYPE basis.  To get the "vectorized" 3x3xN output array, provide 
%   vector inputs in place of scalars, and matrix inputs in place of
%   vectors.  The convention for the latter is that the columns of the
%   input matrices correspond to input vectors for each desired output.
%
%keyword: Coordinate Transformations, Attitude, Utilities,
% See also: JATDCM
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

% Russell Carpenter
% NASA GSFC
% Date: 2005-02-22 12:07:47 -0500 (Tue, 22 Feb 2005)

type = lower(type(1:3));
switch type,
    case {'ax1','ax2','ax3','axs'},
        theta = varargin{1};
        n = length(theta);
        cth = cos(theta(:))';
        sth = sin(theta(:))';
        zero = zeros(1,n);
        one = ones(1,n);
        switch type,
            case 'ax1',
                u1 = [one; zero; zero];
                u2 = [zero;  cth; sth];
                u3 = [zero; -sth; cth];
            case 'ax2',
                u1 = [cth; zero; -sth];
                u2 = [zero; one; zero];
                u3 = [sth; zero;  cth];
            case 'ax3',
                u1 = [ cth; sth; zero];
                u2 = [-sth; cth; zero];
                u3 = [zero; zero; one];
            case 'axs', % Hughes, Spacecraft Attitude Dynamics, p.12
                a = varargin{2};
                mth = one - cth;
                u1 = [ ...
                    mth.*a(1,:).*a(1,:) + cth;
                    mth.*a(1,:).*a(2,:) + a(3,:).*sth;
                    mth.*a(1,:).*a(3,:) - a(2,:).*sth];
                u2 = [ ...
                    mth.*a(1,:).*a(2,:) - a(3,:).*sth;
                    mth.*a(2,:).*a(2,:) + cth;
                    mth.*a(2,:).*a(3,:) + a(1,:).*sth];
                u3 = [ ...
                    mth.*a(1,:).*a(3,:) + a(2,:).*sth;
                    mth.*a(3,:).*a(2,:) - a(1,:).*sth;
                    mth.*a(3,:).*a(3,:) + cth];
        end
    case {'ric','vnb'},
        r = varargin{1};
        v = varargin{2};
        uh = unit(cross(r,v));
        switch type,
            case 'ric',
                u3 = uh;
                u1 = unit(r);
                u2 = cross(u3,u1);
            case 'vnb',
                u2 = uh;
                u1 = unit(v);
                u3 = cross(u1,u2);
        end
    case 'los',
        r1 = varargin{1};
        r2 = varargin{2};
        u1 = unit(r2-r1);
        u2 = unit(cross(r2,r1));
        u3 = cross(u1,u2);
    case 'con',
        dr = varargin{1};
        dv = varargin{2};
        u1 = unit(dr);
        u3 = unit(cross(dr,dv));
        u2 = cross(u3,u1);
end
D(3,:,:) = u3;
D(2,:,:) = u2;
D(1,:,:) = u1;
