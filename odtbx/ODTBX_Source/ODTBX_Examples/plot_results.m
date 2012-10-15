function plot_results(t,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt,S,sca)
% PLOT_RESULTS Plots estimator error, measurment residuals, variance sandpiles 
% 
% The first 5 input parameters are cell arrays of dimension equal to the 
% number of Monte Carlo cases.  In the dimensions below, 
%  k  = no. of time points
%  ns = no. of solve-for states
%  n  = no. of true states
%  m  = no. of measurements
% INPUT:
%  t{i}     : (k x 1) time vector
%  P{i}     : (ns x k) "scrunched" array of estimate error covariance
%  e{i}     : (ns x k) estimate error
%  dy{i}    : (m x k) either measurement innovations or residuals
%  Pa       : (n x n x k) true covariance due to initial error
%  Pv       : (n x n x k) true covariance due to measurement noise
%  Pw       : (n x n x k) true covariance due to process noise
%  Phata    : (ns x ns x k) estimator covariance due to initial error
%  Phatv    : (ns x ns x k) estimator covariance due to measurement noise
%  Phatw    : (ns x ns x k) estimator covariance due to process noise
%  Pdy{i}   : (m x m x k) formal covariance H*P*H'+R corresponding to dy
%  Pdyt     : (m x m x k) true covariance H*P*H'+R corresponding to dy
%  S        : (ns x n) solve-for mapping matrix [default = eye(n)]
%  sca      : scalar value that scales e and the P's [default = 1]
%
%      See also: estval, plot_ominusc, varpiles
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

% Original Creation:
% ---------------------------------
% Sun Hur-Diaz
% Emergent Space Technologies, Inc.
% July 29, 2009
%
% Modification History:
% ---------------------------------
% 2009/09/15  S. Hur-Diaz    Modified to use H*P*H'+R for the measurement
%                            error covariance

mhat=size(Phata);
m=size(Pa);

if nargin < 13 || isempty(S)
    if mhat(1) < m(1)
        error('Need to specify S in the function plot_results')
    else
        S=eye(mhat(1));
    end
end
if nargin < 14
    sca = 1;
end
sca2=sca*sca;

% Scale the errors and error covariances by sca and sca^2 for plotting
for k = 1:length(e),
    P{k} = P{k}*sca2;
    e{k} = e{k}*sca;
end
Pa = Pa*sca2; Pv = Pv*sca2; Pw = Pw*sca2;
Phata = Phata*sca2; Phatv = Phatv*sca2; Phatw = Phatw*sca2;
%
% Compute the true covariance mapped to the solve-for space
for k = length(t{1}):-1:1,
    SPaSt(:,:,k) = S*Pa(:,:,k)*S'; %#ok<AGROW>
    SPvSt(:,:,k) = S*Pv(:,:,k)*S'; %#ok<AGROW>
    SPwSt(:,:,k) = S*Pw(:,:,k)*S'; %#ok<AGROW>
end
SPSt = SPaSt + SPvSt + SPwSt;

% Plot the estimator errors and estimator covariance information
estval(t,e,P,scrunch(Pa+Pv+Pw),0);

% Plot the measurement residuals or innovations
plot_ominusc(t,dy,Pdy,Pdyt,ceil(mhat(1)/3));
%
% Compute the covariance differences for plotting variance sandpiles:
dPa = SPaSt - Phata;
dPv = SPvSt - Phatv;
dPw = SPwSt - Phatw;
Phatlin = Phata + Phatv + Phatw; % Total Phat corresponding to lin cov
for k = 1:length(S),
    figure(gcf+1)
    varpiles(t{1},dPa(k,k,:),dPv(k,k,:),dPw(k,k,:),...
        SPaSt(k,k,:),SPvSt(k,k,:),SPwSt(k,k,:),...
        Phata(k,k,:),Phatv(k,k,:),Phatw(k,k,:),SPSt(k,k,:),Phatlin(k,k,:));
end