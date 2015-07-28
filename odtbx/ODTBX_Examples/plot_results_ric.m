function plot_results_ric(t,x,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,Pdy,Pdyt,S,sca)
% Plots estimator error, measurment residuals, variance sandpiles in the Radial-In-track-Cross-track directions.
%
% PLOT_RESULTS_RIC Plots estimator error, measurment residuals, variance sandpiles 
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
%      See also: estval, plot_results, varpiles
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
    SPaSt(:,:,k) = S*Pa(:,:,k)*S'; 
    SPvSt(:,:,k) = S*Pv(:,:,k)*S'; 
    SPwSt(:,:,k) = S*Pw(:,:,k)*S'; 
end
for k = length(t{1}):-1:1,
    M = dcm('ric',x{1}(1:3,k),x{1}(4:6,k));
    M = blkdiag(M,M);
    SPaSt(1:6,1:6,k) = M*SPaSt(1:6,1:6,k)*M';
    SPvSt(1:6,1:6,k) = M*SPvSt(1:6,1:6,k)*M';
    SPwSt(1:6,1:6,k) = M*SPwSt(1:6,1:6,k)*M';
    Phata(1:6,1:6,k) = M*Phata(1:6,1:6,k)*M';
    Phatv(1:6,1:6,k) = M*Phatv(1:6,1:6,k)*M';
    Phatw(1:6,1:6,k) = M*Phatw(1:6,1:6,k)*M';
end

SPSt = SPaSt + SPvSt + SPwSt;

for el = 1:length(t),
    for k = 1:length(t{el}),
        M = dcm('ric',x{el}(1:3,k),x{el}(4:6,k));
        M = blkdiag(M,M);
        e{el}(1:6,k) = M*e{el}(1:6,k);%*1e3;
        phatel = unscrunch(P{el}(:,k));
        phatel(1:6,1:6) = M*phatel(1:6,1:6)*M';
        phatel = (phatel + phatel')/2;
        P{el}(:,k) = scrunch(phatel);%*1e6;
    end
end
% Plot the estimator errors and estimator covariance information
estval(t,e,P,scrunch(SPSt),0);

igcf = gcf;

for fh = igcf-1:igcf,
    fa = get(fh,'children');
    for k = 1:3,
        axes(fa(k))
        xlabel('Elapsed Time [Days]')
        if fh == igcf-1
            if k == 3
                title('Pos. Error in Radial Dir. [km]','interpreter','none')
            elseif k == 2
                title('Pos. Error in In-Track Dir. [km]','interpreter','none')
            elseif k == 1
                title('Pos. Error in Cross-Track Dir. [km]','interpreter','none')
            end
        else
            if k == 3
                title('Vel. Error in Radial Dir. [km/s]','interpreter','none')
            elseif k == 2
                title('Vel. Error in In-Track Dir. [km/s]','interpreter','none')
            elseif k == 1
                title('Vel. Error in Cross-Track Dir. [km/s]','interpreter','none')
            end
        end
    end
end

% Plot the measurement residuals or innovations
nmeas = size(dy{1},1);
if isempty(Pdy)
    Pdy = deal(cell(1,length(dy)));
    [Pdy{:}] = deal(NaN(nmeas,nmeas,m(3)));
end
if isempty(Pdyt)
    Pdyt = NaN(nmeas,nmeas,m(3));
end
plot_ominusc(t,dy,Pdy,Pdyt,ceil(mhat(1)/3));


% return
%
% Compute the covariance differences for plotting variance sandpiles:
dPa = SPaSt - Phata;
dPv = SPvSt - Phatv;
dPw = SPwSt - Phatw;
Phatlin = Phata + Phatv + Phatw; % Total Phat corresponding to lin cov
    for k = 1:length(S),
    if verLessThan('matlab','8.4.0')
        % execute code for R2014a or earlier
        figure(gcf+1)
    else
        % execute code for R2014b or later
        current_fig = gcf;
        figure(current_fig.Number+1)
    end
    varpiles(t{1},dPa(k,k,:),dPv(k,k,:),dPw(k,k,:),...
        SPaSt(k,k,:),SPvSt(k,k,:),SPwSt(k,k,:),...
        Phata(k,k,:),Phatv(k,k,:),Phatw(k,k,:),SPSt(k,k,:),Phatlin(k,k,:));
end