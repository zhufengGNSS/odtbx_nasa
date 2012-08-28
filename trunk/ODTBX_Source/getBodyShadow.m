function [c_shad] = getBodyShadow(Epoch,satPos,ShadowBodies,aunit)
% getBodyShadow Determine whether any of the specified celestial bodies is 
% shadowing the spacecraft. This function was extracted from srpAccel.m
% in order to make it available to the regression/validation tests.
%
% See also: srpAccel.m
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
%   Ravi Mathur         08/28/2012      Extracted from srpAccel.m

c_shad = 1;

if ~strcmpi(ShadowBodies(1),'NONE')

	for k = 1:length(ShadowBodies)

		Body = char(ShadowBodies(k));

		if strcmpi(Body,'Earth') | strcmpi(Body,'Moon')


			% Determine sun vector relative to body (km)

			rsun_b = ephemDE405('SUN',Epoch,'UTC') - ephemDE405(Body,Epoch,'UTC');

			% Determine satellite position relative to body (km)

			satPos_b = satPos + ephemDE405('EARTH',Epoch,'UTC') - ephemDE405(Body,Epoch,'UTC');

			% planet centric distance of the sun (km)

			rmsun = norm(rsun_b);

			% planet centric unit vector of the sun

			usun = unit(rsun_b);

			% planet centric distance of the spacecraft (km)

			rscm = norm(satPos_b);

			% compute unit position vector of spacecraft

			usat = unit(satPos_b);

			% determine planetary body equitorial radius

			req = JATConstant('meanRadius',Body)/1e3;

			% determine shadow conditions

			a = usat(2) * usun(3) - usat(3) * usun(2);
			b = usat(3) * usun(1) - usat(1) * usun(3);
			c = usat(1) * usun(2) - usat(2) * usun(1);

			d = norm([a  b c]);

			e = dot(usat, usun);

			u = asin(0.00460983743 / (rmsun / aunit));
			p = asin(0.0046951089 / (rmsun / aunit));

			if (e > 0)
				q = -d;
			else
     				q = d;
			end

			if strcmpi(upper(Body),'EARTH')

				% increase the Earth's radius by 90 km
				% to account for the atmosphere

				ratm = req + 90;
			else

				% no atmosphere information is available
				% for other celestial bodies

				ratm = req;

			end

			PenCutoff = ratm/tan(u);
	
			b = asin(ratm / rscm);

			v = b - u;
			w = b + p;

			x = sin(v);
			y = sin(w);
	
			% determine shadow conditions
			% determine c_shad (=0 for in shadow, =1 for in Sunlight)

			if (q <= y & q > x) & c_shad == 1 & rscm < PenCutoff
	     			% penumbra
	     			c_shad = 0.5;
			elseif (q <= x & q >= 0)
	     			% umbra
	     			c_shad = 0;
			else
     				% sunlight
			end

		else

			warning(['Shadow model not available for ',Body])

		end

	end

end

end