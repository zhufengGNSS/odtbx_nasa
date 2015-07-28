function varargout = plot_ominusc(t,dy,Pdy,Pdyt,fhs,iusemn,eflag,sig)
% PLOT_OMINUSC  Plots the observed measurement minus the computed
%
% plot_ominusc(t,dy,Pdy) plots the measurement errors, i.e., the observed
% minus the computed given by the cell dy as dots vs. time given by 
% the cell t.  It also plots the standard deviations of these errors from 
% the measurement error covariance Pdy, which is equal to H*P*H'+R.  These
% are plotted as green plusses.
%
%      The plots are distributed among multiple figure windows, using
% maximum of 3 plots per figure.  If the data represent more than one 
% Monte Carlo case, then it also plots the ensemble mean in a blue line 
% with x's, and the ensemble +/- 1- and 3-sigma deviations of the data are 
% computed and plotted as cyan lines.  For multiple Monte Carlo cases,
% the number of measurement types as well as the time vector have to be the
% same for all the cases.   
%
% NOTE: The measurement error dy can be either innovations or residuals, 
% depending on what is supplied by the user.  It is assumed that Pdy (and 
% Pdyt) correspond to the dy supplied.
%
% plot_ominusc(t,dy,Pdy,Pdyt) also plots true measurement error
% covariance from the true components of H*P*H'+R in black plusses.
%
% plot_ominusc(t,dy,Pdy,Pdyt,fhs) specifies the figure window
% number from which to increment the new figure windows that will be
% created.  Default is to use largest existing figure handle.
%
% plot_ominusc(t,dy,Pdy,Pdyt,fhs,iusemn) specifies whether the plot of the
% ensemble standard deviation of the error is relative to the ensemble mean
% of the error (iusemn=1) or relative to zero (iusemn=0).  Default is
% iusemn=1.
%
% plot_ominusc(t,dy,Pdy,Pdyt,fhs,iusemn,eflag) plots any edited
% measurements as red dots, while keeping the accepted measurements as
% mustard dots. 
%
% plot_ominusc(t,dy,Pdy,Pdyt,fhs,iusemn,eflag,sig) allows the user to set
% what sigma value they wish to plot. Default is set to 2. 
%
% INPUT:
%   t{i}    : m-vector of times for each case i
%   dy{i}   : (n x m) array of measurement erors, n is the no. of 
%             measurements and m is the number of time steps corresponding
%   Pdy{i}  : (n x n x m) array of measurement erors, n is the no. of 
%             measurements and m is the number of time steps corresponding
%   Pdyt{i} : (n x n x m) array of measurement erors, n is the no. of 
%             measurements and m is the number of time steps corresponding
%   fhs     : scalar integer corresponding to the figure window number to
%             increment from
%   iusemn  : scalar indicating whether to plot standard deviations
%             relative to the mean
%   eflag   : (n x m) array of logical elements where 0 (false) represents
%             an edited measurement
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
% Emergent Space Technologies, Inc.
% 2009 July 23
%
% Modification History
% --------------------------
% 2009/09/15  S. Hur-Diaz   Replaced measurement noise R's with measurement
%                           error covariances.  Also, implemented ensemble
%                           mean and std.  Added regression tests.
%
% 2010/08/24  J. Gaebler    Added dynamic legend creation. Added eflag to 
%                           plot edited measurements separately. Inserted 
%                           a warning if the time elements between Monte
%                           Carlo cases does not match. 
%
% 2012/11/04  R. Mathur     Auto-computed figure handle now correctly
%                           ignores existing non-integer figures. This
%                           fixes the incompatibility with 'publish'.

% Self test for regression testing
if nargin < 2, % test data
    my=100; % # of time steps
    nc=10;  % # of cases
    mx=19;  % # of states
    RandStream.setGlobalStream(RandStream('shr3cong', 'Seed', 0))
    e=randn(mx,my,nc); 
%     fs=e+1.5;
    for i=nc:-1:1
        dy{i}=e(:,:,i);
        fs{i}=e(:,:,i)+1.5;
    end
        
    if nargin == 0
        % Case 1
        fail1 = plot_ominusc(0);
        % Case 2
        fail2 = plot_ominusc(1);
        varargout{1} = any([fail1 fail2]);
        return
    else
        iusemn = t;
        icase = iusemn;
        t={linspace(0,100,my)};
    end
elseif nargin >= 2,
    nc = length(dy);       % number of Monte Carlo cases
    [mx,my] = size(dy{1});  % mx = no. of states, my = number of time steps
    
    % Determine the sigmas from Pdy
    fs=cell(1,nc);
    if nargin >=3 && ~isempty(Pdy)
        for i=1:nc
            for j=my:-1:1
                fs{i}(:,j)=sqrt(diag(Pdy{i}(:,:,j)));
            end
        end
    else
        fs = cell(1,nc);
        for i=1:nc
            fs{i} = NaN(size(dy{i}));
        end
    end

    if nargin<6 || isempty(iusemn)
        iusemn=1;
    end
    
    if nargin<7 || isempty(eflag)
        e_ind=cell(nc,1);
        for i=1:nc
            e_ind{i}=(true(size(dy{i})));
        end
    else
        e_ind=eflag;
    end
    
    if nargin<8 || isempty(sig)
        sig=2;
    end
end

nf=3;                   % Number of plots in one figure
nsub=min(nf,mx);        % Maximum no. of plots in a figure window

% True sigmas from the true measurement error covariance
if nargin<4 || isempty(Pdyt)
    ts = NaN(mx,my);
else
    for j=my:-1:1
        ts(:,j) = sqrt(diag(Pdyt(:,:,j)));
    end
end

nfig = ceil(mx/nf);      % Number of figure windows needed
marksize = 5;
linewidth = 2;

color_data = [.8 .8 .2];

if nc > 1
    if iusemn == 1
        str_estd = 'Ensemble 1- to 3-STD of dy wrt ensemble mean';
    elseif iusemn == 0
        str_estd = 'Ensemble 1- to 3-STD of dy wrt zero';
    else
        error('Need to specify 1 or 0 for the iusemn input parameter in plot_ominusc.m')
    end
    
    stat_flag=true;
    try
        e = reshape([dy{:}],mx,my,nc);
        t = [t{1}];          % same time steps for all cases
        fs=reshape(cell2mat(fs),mx,my,nc);
    catch
        warning('Each Monte Carlo case must have the same time steps to calculate ensemble statistics.')
        stat_flag=false;
    end
end

% If a starting figure handle was not specified, then set it to the
% largest existing figure handle.
if(~exist('fhs', 'var') || isempty(fhs))
    % There may be figures with non-integer handles, e.g. if the
    % 'publish' command was used. We need to eliminate these.
    if verLessThan('matlab','8.4.0')
        % execute code for R2014a or earlier
        allfigs = [get(0, 'children');0]; % Get all open figures
        badidx = find(allfigs ~= floor(allfigs)); % Find non-integer handles
        allfigs(badidx) = []; % Remove non-integer handles
        fhs = max(allfigs); % Get largest of remaining integer handles 
    else
        % execute code for R2014b or later
        allfigs = [get(groot, 'children')];% Get all open figures
        allfigs = allfigs.Number;%Convert object to number
        badidx = find(allfigs ~= floor(allfigs)); % Find non-integer handles
        allfigs(badidx) = []; % Remove non-integer handles
        fhs = max(allfigs); % Get largest of remaining integer handles 
    end
    
end

if exist('fhs') 
    if verLessThan('matlab','8.4.0')
        % execute code for R2014a or earlier

    else
        % execute code for R2014b or later
%        isfhs = isnumeric(fhs)
        if isnumeric(fhs) == 0;
            fhs = fhs.Number;
        end
    end
end

for j=1:nfig
    ifig=j+fhs;    % Determine figure number
    figure(ifig)

    for i=1:3
        ip=(j-1)*nf+i;
        if ip>mx, break, end
        subplot(nsub,1,i)

        if nc > 1 % More than 1 Monte Carlo case
            if stat_flag
                em = NaN(mx,my);        % Ensemble mean of dy
                es = NaN(mx,my);        % Ensemble std of dy
                fsem = NaN(mx,my);      % Ensemble mean of sigmas from the formal Pdy
                tr = repmat(t(:)',3,1);
                
                % Find the mean and the standard dev of non-NaN data
                for k=1:my
                    q = ~isnan(e(ip,k,:));
                    em(ip,k) = mean(e(ip,k,q));
                    es(ip,k) = std(e(ip,k,q));
                    fsem(ip,k) = mean(fs(ip,k,q));
                end
                
                % Form the 1-sig and the 3-sig data for plot
                esp  = reshape([iusemn*[1;1]*em(ip,:)+[1;3]*es(ip,:); NaN(1,my)],3*my,1);
                esm  = reshape([iusemn*[1;1]*em(ip,:)-[1;3]*es(ip,:); NaN(1,my)],3*my,1);
                
                hh1 = plot(t,squeeze(e(ip,:,:)),'y.','color',color_data);  % meas dy
                hold on
               
                hh2 = plot(t,em(ip,:),'bx-','linewidth',linewidth);        % ensemble mean dy
                
                hh3a = plot(tr(:),esp,'c','linewidth',linewidth);          % + sig of ensemble dy
                hh3b = plot(tr(:),esm,'c','linewidth',linewidth);          % - sig of ensemble dy
                
                hh4a = plot(...                                            % STD dy
                    t,2*fsem(ip,:),'g+',...                                % 2sig formal
                    t,2*ts(ip,:),'k+');                                    % 2sig true
                hh4b = plot(...                                            % -STD dy
                    t,-2*fsem(ip,:),'g+',...                               % -2sig formal
                    t,-2*ts(ip,:),'k+');                                   % -2sig true
                
                hold off
                h=[hh1;hh2;hh3a;hh3b;hh4a;hh4b];
                
                if i==1 %(t,dy,Pdy,Pdyt,fhs,iusemn,eflag)
                    if nargin>=4 && ~isempty(Pdyt)
                        legend([hh1(1);hh2;hh3a;hh4a],'Measurement Error dy','Ensemble Mean of dy',str_estd,...
                            '2X ensemble mean of sigmas from formal HPH''+R','2-sigma from true HPH''+R')
                    elseif ~all(isnan(fs))
                        legend([hh1(1);hh2;hh3a;hh4a(1)],'Measurement Error dy','Ensemble Mean of dy',str_estd,...
                            '2X ensemble mean of sigmas from formal HPH''+R')
                    else
                        legend([hh1(1);hh2;hh3a],'Measurement Error dy','Ensemble Mean of dy',str_estd)
                    end
                    title('Observed Minus Computed Measurements')
                end
            else % this is where plotting of MC cases with dissimilar time steps will go
            end
            
        else % Only one case. Cell has one element
            h1 = plot(t{1}(e_ind{1}(ip,:)~=0),dy{1}(ip,e_ind{1}(ip,:)~=0),'y.',...
                t{1},squeeze(sig*fs{1}(ip,:,1)),'g+',...  % sig formal
                t{1},sig*ts(ip,:),'k+',...             % sig true
                t{1}(e_ind{:}(ip,:)==0),dy{1}(ip,e_ind{:}(ip,:)==0),'r.'); % edited measurments
            hold on
            h2 = plot(t{1},-squeeze(sig*fs{1}(ip,:,1)),'g+',... % -sig formal
                t{1},-sig*ts(ip,:),'k+');              % -sig true
            h=[h1;h2];
            set(h(1),'color',color_data) % change to mustard color
            if i==1
                leg_h=h1;
                leg_str={'Measurement Error dy',[num2str(sig),'-sigma from formal HPH''+R'],...
                    [num2str(sig),'-sigma from true HPH''+R'],'Measurements edited'};
                if nargin<7 || isempty(eflag),leg_str(4)=[];end
                if nargin<4 || isempty(Pdyt),leg_str(3)=[];leg_h(3)=[];end
                if nargin<3 || isempty(Pdy),leg_str(2)=[];leg_h(2)=[];end
                legend(leg_h,leg_str)
                title('Observed Minus Computed Measurements')
            end
        end
                
        set(h,'markersize',marksize)
        set(gca,'box','off','fontname','times','fontsize',11,...
            'xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1])
        zoom yon
        if ip-(j-1)*nf~=nsub,
            set(gca,'xtick',[],'xcolor',get(gca,'color'))
        end
        ylabel(['dy_{',num2str(ip),'}'],'fontname','times','fontsize',12,...
            'fontangle','italic','rotation',0,...
            'horizontalalign','right','verticalalign','middle')
    end
end
for i = fhs+ (1:nfig),
    set(i,'color',get(gca,'color'),'paperposition',[.25 .5 8 10])
end

% Regression testing
if nargin < 2 && nargout > 0
    load plot_ominusc_test_data
    lh = length(h);
    fail = 0;
    if icase==0
        g_save = g_save0;
    else
        g_save = g_save1;
    end
    for i = lh:-1:1
        g(i) = get(h(i));
        fail = any([fail ( any (abs(g(i).YData - g_save(i).YData) > eps) )]);
    end
    varargout{1} = fail;
    eval(['close ',num2str(fhs+(1:nfig))]);
end