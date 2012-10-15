function Sta = stationData(ID,Sta)

% STATIONDATA Provides database information about selected ground stations
%
%  Sta = stationData(ID,Sta)
%
%  This function returns information about the selected ground stations from
%  table data.  The data includes ground station position information and
%  average monthly climate information.
%
%  Note that the average monthly climate information is specific to three
%  ground stations (12, 41, & 61).  These three sites are located in California,
%  Australia, and Spain respectively.  Since all of the ground stations are
%  located in these three regions, the climate data for these three ground
%  stations are used for all of the stations in that region.
%
%  There are three ways to call the staionData function.
%
%  1. Sta = stationData;        - Provides all station data
%  2. Sta = stationData(ID)     - Provides data for IDs specified in ID vector
%  3. Sta = stationData(ID,Sta) - Adds to existing Sta structure
%                                   (overwriting only the fields provided)
%
% INPUTS
%      VARIABLE         SIZE        DESCRIPTION (Optional/Default)
%      ID               (1xN)       (Optional) Vector of desired station IDs
%      Sta              Structure   (Optional) Structure to append with 
%                                     provided fields
%
% OUTPUTS
%      Sta.DSS_IDs      (1xN)       Station IDs (matches ID input)   
%      Sta.siteIDs      (1xN)       Site locations (1=California,
%                                     2=Australia, 3=Spain)
%      Sta.lat          (1xN)       Station Latitude (rad)
%      Sta.lon          (1xN)       Station Longitude (rad)
%      Sta.height       (1xN)       Station Height (m)
%      Sta.staPos       (3xN)       Station ECEF Position (m)
%      Sta.posUnc       (3xN)       Station Position Uncertainty (m)
%      Sta.Press        (1xN)       Average monthly pressure (bar)
%      Sta.Temp         (1xN)       Average monthly temperature (K)
%      Sta.RelHum       (1xN)       Average monthly relative humidity (%)
%      Sta.TempLapse    (1xN)       Average monthly temperature lapse (K/m)
%                                     (-delta Temp / delta altitude)
%
% Available station IDs
%
% [12 13 14 15 16 17 23 24 25 26 27 28 33 34 42 43 45 46 53 54 61 63 65 66]
%
% NOTES
%
%   This code is based on Will Campbell's work.
%
% VALIDATION/REGRESSION TEST
%
%   These tests have been moved to stationData_test.m to conform to the
%   new regression testing format.
%
%   keywords: ground station
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

%  REVISION HISTORY
%   Author      		Date         	Comment
%               	   (MM/DD/YYYY)
%   Keith Speckman      06/10/2008          Original
%   Ravi Mathur         08/28/2012      Extracted regression test

% Data monthly averages for DSN Stations 12, 41, & 61

% STATION LOCATIONS
%   "stn"  - vector entry number of antenna
%   "site" - vector entry number of site
% Station IDs (each antenna's DSS Number)
DSS_IDs = [12 13 14 15 16 17 23 24 25 26 27 28 33 34 42 43 45 46 53 54 61 63 65 66];
% Site IDs (1=California 2=Australia 3=Spain)
site_IDs = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3];
% STATION LAT-LON-ALT (rad-rad-km)
gs_LLA = [0.61610013885628   4.24454497545834   0.962875;
   0.61517907843179   4.24473748508091   1.071178;
   0.61829862059561   4.24307803931200   1.002114;
   0.61822797703666   4.24311893470348   0.973945;
   0.61682623528894   4.24335534575628   0.944711;
   0.61686644510540   4.24335786581780   0.937650;
   0.61679151991442   4.24336907664936   0.946086;
   0.61679749378860   4.24333537235377   0.952156;
   0.61675768549484   4.24332543338241   0.960862;
   0.61672412489802   4.24336638661217   0.970159;
   0.61502387551287   4.24504830281645   1.053203;
   0.61502387958530   4.24500919493301   1.065382;
  -0.60387538177882   2.60024544079269   0.684839;
  -0.60391039225943   2.60022576298108   0.692750;
  -0.60387206778642   2.60021359784226   0.675356;
  -0.60384152966354   2.60021359774530   0.689608;
  -0.60391076047542   2.60015108577230   0.675086;
  -0.60379638992221   2.60024526465988   0.677551;
   0.70559046406420   6.20901484054325   0.827501;
   0.70553846731208   6.20895548155602   0.823939;
   0.70561458218735   6.20902583640843   0.841159;
   0.70565770151615   6.20904352808384   0.865544;
   0.70558745914052   6.20898402166477   0.834539;
   0.70563614956029   6.20898402733709   0.850582];
% STATION ECEF XYZ POSITION (km)
gs_ECEF = [-2350.443812 -4651.980837  3665.630988;
    -2351.112491 -4655.530714 3660.912787;
    -2353.621251 -4641.341542 3677.052370;
    -2353.538790 -4641.649507 3676.670043;
    -2354.763158 -4646.787462 3669.387069;
    -2354.730357 -4646.751776 3669.440659;
    -2354.757567 -4646.934675 3669.207824;
    -2354.906495 -4646.840128 3669.242317;
    -2355.022066 -4646.953636 3667.040895;
    -2354.890967 -4647.166925 3668.872212;
    -2349.915260 -4656.756484 3660.096529;
    -2350.101849 -4656.673447 3660.103577;
    -4461.083514 2682.281745 -3674.570392;
    -4461.146756 2682.439293 -3674.393542;
    -4460.981016 2682.413525 -3674.582072;
    -4460.894585 2682.367554 -3674.748580;
    -4460.935250 2682.765710 -3674.381402;
    -4460.828619 2682.129556 -3674.975508;
     4849.330129 -0960.338092 4114.758766;
     4849.519988 -0360.641653 4114.504590;
     4849.245211 -0360.278166 4114.884445;
     4849.092647 -0360.180569 4115.109113;
     4849.336730 -0360.488859 4114.748775;
     4849.148543 -0360.474842 4114.995021];
% STATION ECEF XYZ POSITION UNCERTAINTY (km)
gs_unc = [0.0225    0.0446    0.0366;
    0.0300    0.0594    0.0476;
    0.0187    0.0368    0.0304;
    0.0179    0.0352    0.0287;
    0.0448    0.0883    0.0706;
    0.0448    0.0883    0.0706;
   20.7926   41.0326    0.0699;
    0.0448    0.0883    0.0706;
    0.0448    0.0883    0.0706;
    0.0448    0.0883    0.0706;
    0.0446    0.0884    0.0706;
    0.0446    0.0884    0.0706;
    0.0851    0.0512    0.0699;
    0.0860    0.0517    0.0706;
    0.0485    0.0292    0.0387;
    0.0394    0.0237    0.0310;
    0.0366    0.0220    0.0289;
    0.0439    0.0264    0.0351;
    0.0950    0.0071    0.0761;
    4.4801    0.3332    3.5822;
    0.0499    0.0037    0.0380;
    0.0440    0.0033    0.0327;
    0.0424    0.0032    0.0310;
    0.0950    0.0071    0.0761]/1E3;

% Atmospheric Data for DSN Ground Stations
%
% Data monthly averages for DSN Stations 12, 41, & 61
% SHOULD ADJUST TO REFLECT ACTUAL TRACKING STATIONS

P = [906.72 992.46 931.37;
    906.66 0994.72 928.26;
    902.70 1003.85 926.73;
    900.03 1001.52 923.12;
    898.66 1003.90 924.43;
    895.46 1005.05 927.58;
    900.91 1000.80 926.55;
    900.54 0936.61 925.93;
    899.64 1001.92 926.18;
    904.38 1000.70 927.74;
    903.62 0997.60 925.67;
    902.58 0996.81 929.24]/1000;

T = [20.08 34.94 12.14;
    18.12 26.68 12.14;
    16.18 21.44 16.05;
    16.16 28.66 13.97;
    20.78 18.92 15.09;
    22.96 20.13 20.08;
    30.50 13.70 27.13;
    31.60 14.00 24.13;
    23.41 18.68 21.05;
    27.96 21.37 19.07;
    19.28 19.70 13.92;
    12.54 24.97 11.42] + 273;

RH = [.1930 .2839 .5267;
    .2783 .3029 .6245;
    .2748 .4555 .5092;
    .2504 .3621 .6705;
    .2978 .4709 .6073;
    .3070 .4201 .5280;
    .3312 .5316 .3927;
    .3561 .4622 .4857;
    .4384 .3418 .5718;
    .2129 .3137 .5738;
    .2814 .3662 .7404;
    .3026 .3669 .6041];

GAM = [6.20 7.26 6.66;
    6.74 5.79 6.60;
    6.53 5.59 6.53;
    6.08 7.10 6.72;
    7.01 5.86 5.84;
    6.40 6.28 6.25;
    7.11 5.50 6.55;
    7.53 5.21 6.23;
    6.25 5.61 6.01;
    6.23 5.75 5.95;
    6.00 5.36 6.04;
    6.18 6.14 6.09];

if nargin == 0
% Output all station data

	Sta.DSS_IDs = DSS_IDs;
	Sta.siteIDs = site_IDs;
	Sta.lat = gs_LLA(:,1)';
	Sta.lon = gs_LLA(:,2)';
	Sta.height = gs_LLA(:,3)'*1e3;
	Sta.staPos = gs_ECEF'*1e3;
	Sta.posUnc = gs_unc'*1e3;
	Sta.Press = P(:,site_IDs);
	Sta.Temp = T(:,site_IDs);
	Sta.RelHum = RH(:,site_IDs);
	Sta.TempLapse = GAM(:,site_IDs)/1e3;

else
% Output data for specified IDs

	for k = 1:length(ID)
		StaIndexTemp = find(DSS_IDs==ID(k));
		if isempty(StaIndexTemp)
			error(['DSS ID ',num2str(ID(k)),' is not available'])
		end
		StaIndex(k) = StaIndexTemp;
	end

	Sta.DSS_IDs = ID;
	Sta.siteIDs = site_IDs(StaIndex);
	Sta.lat = gs_LLA(StaIndex,1)';
	Sta.lon = gs_LLA(StaIndex,2)';
	Sta.height = gs_LLA(StaIndex,3)'*1e3;
	Sta.staPos = gs_ECEF(StaIndex,:)'*1e3;
	Sta.posUnc = gs_unc(StaIndex,:)'*1e3;
	Sta.Press = P(:,site_IDs(StaIndex));
	Sta.Temp = T(:,site_IDs(StaIndex));
	Sta.RelHum = RH(:,site_IDs(StaIndex));
	Sta.TempLapse = GAM(:,site_IDs(StaIndex))/1e3;

end

end

