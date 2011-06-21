function [y,z]= testeom(t,x)

%TESTEOM Example EOM file for use with the MATLAB Adaptor.
%   TESTEOM example EOM file that gives an example of how to input a 
%   user-defined EOM file into the Matlab Adaptor.
%  
%   INPUTS
%   VARIABLE    SIZE    DESCRIPTION (Optional/Default)
%      t        (1xN)     Input to function, time in seconds
%      x        (1xN)      '', state
%
%   OUTPUTS 
%      y        (1xN)     Ouput from function, new state
%      z        (1x1)      '', optional output                                     
%
%    keyword: JAT Adaptor, Example Programs
%    See also jatRK8
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

%   REVISION HISTORY
%   Author      		Date         	Comment
%               		(MM/DD/YYYY)
%   Dave Gaylor     	2006   		Original
%   Kathryn Bradley	04/17/2006		Ammended EOM example & documented
%   Emergent Space Technologies
%   Allen Brown         02/11/2009  updated comments

[m,n] = size(x);
%if m==1 & n==1,
%    warning('DERIVS:scalar', 'Input is a scalar.')
%end

% if x is a cell array, convert to matrix
if iscell(x),
   x = cell2mat(x);
end

% compute eom
y(1) = x(2);
y(2) = -x(1);

z = 0.333;