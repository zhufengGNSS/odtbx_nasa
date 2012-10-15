function [fo,fdo,tab] = ddhermite(xo,x,f,fd,ord)

% DDHERMITE interpolates based on modified divided difference.
%
% [fo,fdo,tab] = ddhermite(xo,x,f,fd,ord)
%   Forms the Hermite interpolating polynominal by modifying the divided 
%   difference basis of order ord and interpolates for the given xo.  The
%   method is based on the course notes by Dr. David Hill of Temple
%   University.
%   (http://astro.temple.edu/~dhill001/course/numanalspring2010/Hermite%20I
%   nterpolation%20Section%205_7.pdf)
%
%   This function was originially created as a prototype of what's
%   implemented in JAT for use by the jatworldpropagator functions.  
%   xo is scalar in this implementation.
% NOTE ddHermite requires all arrays be the same length and equal to ord
% Author Sun Hur-Diaz, Emergent Space Technologies
% Date   02/09/2010
%
%   INPUTS 
%   VARIABLE    SIZE    DESCRIPTION (Optional/Default)
%       xo      (scalar)    currentTime
%       x       (N*1)       time array
%       f       (N*1)       position array
%       fd      (N*1)       velocity array
%       ord     (scalar)    order
%
%   OUTPUTS 
%       fo      (scalar)    position
%       fdo     (scalar)    velocity
%       tab     (ord*2,ord*2+1) divided difference table
%
%   keyword: integrator
%   See also jatWorldPropagatorRK4, jatWorldPropagatorRK8, JATFORCES, 
%   CREATEJATWORLD
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
%               	   (MM/DD/YYYY)
%   Sun Hur-Diaz        02/09/2010   	Original
%   Stephen D Metcalfe  06/09/2010   	Documented and added test for array
%                                       size to be equal to ord in order to
%                                       satisfy requirements of reshape.

lx =length(x);
if length(f) ~= lx || length(fd) ~= lx
    error('Lengths of x, f and fd need to be same.');
end

if nargin < 5 || isempty(ord)
	ord = lx;
elseif ord ~= lx
    error('Lengths of x, f and fd do not match ord.\n');
end

% Form the divided difference table
tab = zeros(2*ord,2*ord+1);
tab(:,1:2) = [reshape([x(:) x(:)]',2*ord,1) reshape([f(:) f(:)]',2*ord,1)]; 
tab(1:2:2*ord-1,3) = fd(:);
tab(2:2:2*(ord-1),3) = (tab(3:2:2*ord-1,2)-tab(2:2:2*ord-2,2))./(tab(3:2:2*ord-1,1)-tab(2:2:2*ord-2,1));
for j = 4:2*ord+1
    di = j-2;
    for i = 1:2*(ord+1)-j
        tab(i,j) = (tab(i+1,j-1)-tab(i,j-1))/(tab(i+di,1)-tab(i,1));
    end
end

% Interpolate on the function
fo = tab(1,2);
p = 1;
v = zeros(2*ord,1);
for i = 3:2*ord+1
    v(i-2) = (xo - tab(i-2,1)); 
    p = p * v(i-2);
    fo = fo + tab(1,i)*p;
end

% Interpolate on the function derivative
fdo = tab(1,3);
for i = 4:2*ord+1
    s = 0;
    for j = 1:i-2
        p = 1;
        for k = [1:j-1 j+1:i-2]
            p = p * v(k);
        end
        s = s + p; 
    end
    fdo = fdo + tab(1,i) * s;
end

