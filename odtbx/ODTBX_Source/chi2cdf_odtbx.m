function [p,dpdx] = chi2cdf_odtbx(x,v)
% ODTBX version of chi2cdf
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
% July 19, 2010

lx = length(x);
p = NaN(size(x));
for i=1:lx
    p(i) = quadgk(@(x)chi2pdf_odtbx(x,v),0,x(i));
end
if nargout>1
    dpdx=chi2pdf_odtbx(x,v);
end