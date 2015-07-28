function h = estval(t,e,P,Pt,fhs,iusemn)
% ESTVAL  Evaluation of estimator outputs (plotting).
%   ESTVAL(T,E,P) plots the estimation errors E (blue lines) and their 
%   "formal" 2-sigma envelope (green dashed lines), derived from the 
%   estimator covariance matrix P, vs. time, T.  ESTVAL automatically 
%   distributes the plots among multiple figure windows, using a maximum 
%   of 3 plots per figure.  
%      If the inputs represent multiple Monte Carlo cases, ESTVAL plots the
%   estimation errors as yellow dots and their ensemble means as blue 
%   lines.  It also plots the ensemble means of their formal sigmas times 2
%   to represent the 2-sigma values (green plusses), as well as a 
%   representation of the actual ensemble deviates of the estimation 
%   errors for each time sample.  This representation is a pair of vertical 
%   cyan lines, with the endpoints of the upper line representing the +1 
%   and +3 sigma actual ensemble standard deviations, and the endpoints of 
%   the lower line the -1 and -3 sigma values:
%
%      ||  <- Top of this line = +3 sigma actual ensemble deviation 
%      ++  <- Positive ensemble mean of "formal" sigmas from P (green plus)
%      ||  <- Bottom of this line = +1 sigma actual ensemble deviation
%      --  <- Ensemble mean of estimation errors E (blue line)
%      ||  <- Top of this line = -3 sigma actual ensemble deviation
%      ++  <- Negative ensemble mean of "formal" sigmas from P (green plus)
%      ||  <- Bottom of this line = -1 sigma actual ensemble deviation
%
%
%   ESTVAL(T,E,P,Pt) also plots the standard deviation derived from the
%   true covariance matrix Pt as black plusses.
% 
%   ESTVAL(T,E,P,Pt,fhs) specifies the figure window number from which to 
%   increment the new figure windows that will be created.  Default is
%   zero.
%
%   ESTVAL(T,E,P,Pt,fhs,iusemn) specifies whether the plot of the ensemble
%   standard deviation of the error is relative to the ensemble mean of the
%   error (iusemn=1) or relative to zero (iusemn=0).  Default is iusemn=1.
%
%   h = ESTVAL(T,E,fs,Pt,fhs,iusemn) assumes that the covariance standard
%   deviations are known directly and given in fs. In this case, fs must be
%   be a 3D matrix. The generated figure handles are returned in h. This
%   mode exists mainly for regression testing.
%
%   N.B.: the current version requires that each M.C. case have the same
%   number of elements, i.e. the number of dimensions and the number of
%   time samples cannot change between cases.  Future versions are
%   anticipated to support differing sample times with an interpolation and
%   extrapolation capability.
%
%   keyword: Estimation, 
%   See also
%      ODEAS estimators:      ESTBAT, ESTSEQ
%      covariance storage:    SCRUNCH, UNSCRUNCH, SIGEXTRACT
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
% NASA Goddard Space Flight Center

% Modification History:
% ---------------------
% 2009 July 21  S. Hur-Diaz
% Added a legend to the first plot, reordered the data plotted so that the 
% lines can be seen better over the error data points, and made the true
% covariance and figure number inputs to be optional.
%
% 2009 Sept 14  S. Hur-Diaz
% Added the option to plot the ensemble sigmas relative to the ensemble
% mean or zero.  Also added a regression test.
%                                                                                                                                                                    -
% 2010 Mar 10 John Gaebler                                                                                                                                           -
%   Added label to plot legend if statistics toolbox and Monte Carlo cases                                                                                           -
%   are present.
%
% 2013 May 3 Ravi Mathur
% Fully extracted internal self-tests to estval_test.m to conform to the
% new regression testing framework.
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
markersize = 5; % Desired size of markers

if nargin < 3, % Old selftest, no longer supported
    error('estval no longer supports direct regression testing. Please use estval_test.');
else
    if iscell(e),
        K = length(e);       % number of Monte Carlo cases
        [n,N] = size(e{1});  % n = no. of states, N = number of time steps
        t = [t{1}];          % same time steps for all cases
        e = reshape([e{:}],n,N,K);
    else
        [n,N,K] = size(e);
    end
    if iscell(P),
        P = reshape([P{:}],n*(n+1)/2,N,K);
        fs = sigextract(P);
    else
        fs = P; % Assume passed-in P was actually standard deviations
    end
    if nargin<6
        iusemn=1;
    end
end
if nargin<4 || isempty(Pt)
    ts = NaN(n,N);
else
    ts = sigextract(Pt);
end
if nargin<5 % && ~exist('fhs','var')
    fhs = 0;
end

if iusemn == 1
    str_estd = 'Ensemble 1- to 3-std of e wrt ensemble mean';
elseif iusemn == 0
    str_estd = 'Ensemble 1- to 3-std of e wrt zero';
else    
    error('Need to specify 1 or 0 for the iusemn input parameter in estval.m')
end

t = t(:)';
em = mean(e,3);          % Ensemble mean of the errors
fsem = mean(fs,3);       % Ensemble mean of the formal sigmas
es = std(e,0,3);         % Ensemble STD of the errors
tr = repmat(t(:)',3,1);  
nax = 3;                 % Maximum number of plots per figure
nfig = ceil(n/nax);      % Number of figure windows needed
% scrn = get(0,'monitorpositions');
% scrn = scrn(end,:);
% fsiz = [min(scrn(3)./[nfig 2]) scrn(4)/2];

for i = 1:n, % Cycle through each state
    ifig = ceil(i/nax);
    figure(ifig+fhs)
%     if strcmp(get(0,'DefaultFigureWindowStyle'),'normal'),
%         set(gcf,'position',[(ifig-1)/nfig*scrn(3) scrn(4)/4 fsiz]);
%     end
    if ifig == nfig,
        iax = rem(n,nax);
    else
        iax = nax;
    end
    if iax == 0,
        iax = nax;
    end
    
    subplot(iax,1,i-(ifig-1)*nax)
    
    if K > 1, % More than one Monte Carlo case
        % Form the 1-sig and 3-sig ensemble STDs
        les  = reshape([iusemn*[1;1]*em(i,:)+[1;3]*es(i,:); NaN(1,N)],3*N,1);
        lesm = reshape([iusemn*[1;1]*em(i,:)-[1;3]*es(i,:); NaN(1,N)],3*N,1);
        h = plot(...
            t,squeeze(e(i,:,1)),'y-',...     % Actual errors
            t(1),em(i,1),'b-',...            % Ensemble mean of e
            tr(1),les(1),'c-',...            % Ensemble 1,3sig of e
            t(1),2*fsem(i,1),'g+-',...        % Ensemble 2sig from P
            t(1),2*ts(i,1),'k+-',...          % 2sig from P true
            tr(:),les,'c-',...               % Ensemble 1,3sig of e
            tr(:),lesm,'c-',...              % Ensemble -(1,3sig) of e
            t,squeeze(e(i,:,2:K)),'y.-',...   % Rest of actual errors
            t,em(i,:),'b-',...               % Ensemble mean of e
            t,2*fsem(i,:),'g+-',...           % Ensemble 2sig from P
            t,2*ts(i,:),'k+-',...             % 2sig from P true
            t,-2*fsem(i,:),'g+-',...          % Ensemble -2sig from P
            t,-2*ts(i,:),'k+-'...            % -2sig from P true
            );
        
        set(h(strcmp(get(h,'marker'),'.')),'color',[.8 .8 .2])
        
        leghan=[h(1:5)];
        legstr={'Actual errors e';'Ensemble mean of e';...
            str_estd;'2X ensemble mean of sigmas from P formal';...
            '2-sigma from P true'};
        
        if ~isempty(ver('stats')) && license('checkout','statistics_toolbox'); % Checks if statistics toolbox is present
            % The above test may not be the best method. Licensing of the 
            % toolboxes may adversely affect the test.
        for j=N:-1:1,
            [~,~,ci(i,j,:)] = vartest(squeeze(e(i,j,:)),ts(i,j)^2);
        end
        hold on
        v = plot(t,[diag([2;2]);diag([-2;-2])]*squeeze(sqrt(ci(i,:,:)))','r+-');
        line(repmat(t,2,1),diag([2;2])*squeeze(sqrt(ci(i,:,:)))','color',[1 .5 .5])
        line(repmat(t,2,1),diag(-[2;2])*squeeze(sqrt(ci(i,:,:)))','color',[1 .5 .5])
        set(gca,'children',circshift(get(gca,'children'),length(h)))
        hold off
        leghan(end+1)=v(1);
        legstr{end+1}='2X confidence interval on ensemble sigmas';
        set(v,'markersize',markersize)
        end
        if i==1 && i-(ifig-1)*nax == 1
            if nargin<4 || isempty(Pt)
                leghan(5)=[];
                legstr(5)=[];
            end
            legend(leghan,legstr);
        end
        
    else % Just one case
        h = plot(t,em(i,:),'b-',...
            t,2*fs(i,:),'g--',...
            t,2*ts(i,:),'k--',...
            t,-2*fs(i,:),'g--',...
            t,-2*ts(i,:),'k--');
        if i-(ifig-1)*nax == 1
            if nargin>=4 && ~isempty(Pt)
                legend('Actual errors','2-sigma from P formal',...
                    '2-sigma from P true');
            else
                legend('Actual errors','2-sigma from P formal');
            end
        end
    end
    set(h,'markersize',markersize)
    ylabel(['x_{',num2str(i),'}'],'fontname','times','fontsize',12,...
        'fontangle','italic','rotation',0,...
        'horizontalalign','right','verticalalign','middle')
    set(gca,'box','off','fontname','times','fontsize',11,...
        'xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1])
    zoom yon
    if i-(ifig-1)*nax~=iax,
        set(gca,'xtick',[],'xcolor',get(gca,'color'))
    end
end
for i = fhs+ (1:nfig),
    set(i,'color',get(gca,'color'),'paperposition',[.25 .5 8 10])
end

kds=get(fhs+1,'children');
set(fhs+1,'children',circshift(kds,length(kds)-find(strcmp(get(kds,'tag'),'legend'))+1));

end % function