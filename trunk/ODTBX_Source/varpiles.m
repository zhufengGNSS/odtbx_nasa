function varargout = varpiles(varargin)
% VARPILES  Variance Sandpiles.
%
%   varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat) plots a
%   stacked area chart showing the time series of the contributions to the
%   total variance of a particular solve-for variable due to _a priori_
%   uncertainty, measurement noise, and process noise.  Both formal and
%   true variance contributions are plotted.  The way in which they are
%   plotted depends on their values, as described below.  The inputs are
%   defined as follows:
%      t     - Time-tags of the variance arrays
%      dVa   - Difference between true and formal a priori variance
%      dVv   - Difference between true and formal measurement noise variance
%      dVw   - Difference between true and formal process noise variance
%      Va    - True a priori variance
%      Vv    - True measurement noise variance
%      Vw    - True process noise variance
%      Vhata - Formal a priori variance
%      Vhatv - Formal measurement noise variance
%      Vhatw - Formal process noise variance
%      V     - Total true variance
%      Vhat  - Total formal variance
%      (OPTIONAL)
%      dVm   - Difference between true and formal maneuver variance
%      Vm    - True maneuver execution variance
%      Vhatm - Formal maneuver execution variance
%
% There are several ways we can plot the sandpiles. When all the delta
% variances are positive, the sandpile will show the formal variance, and
% the deltas due to each component.  When all the deltas are negative, the
% sandpile will show the true variance, and negatives of the deltas.
% Otherwise, varpiles will plot the components of the true variance as a
% positive sandpile, and the components of the formal variance as a
% negative sandpile, and relabel the negative y-axes to indicate this.
%
%   Examples
%
%   Seach Categories (keywords): Estimation
%
%   See also
%      evaluating solutions:  ESTVAL
%      covariance storage:    SCRUNCH, UNSCRUNCH
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

% Modification History
% ---------------------
% 2010/10/18 K Getzandanner Added optional inputs for maneuver execution
%                           variance partitions
% 
% 2011/01/10 K Getzandanner Added logic checks for empty legend entries

%% VARPILES: Variance Sandpiles
%
% "Variance sandpiles" are stacked area charts showing the the time series
% of each solve-for variance's contribution from _a priori_ error variance,
% measurement noise variance, and process noise variance.  
%
% The following mathematical specifications were published from comments
% embedded within the m-file.
%
%% Specifications
%
% The true variance and the formal variance are
%
% $$ P = P_a + P_v + P_w, \quad \hat{P} = \hat{P}_a + \hat{P}_v $$
%
% and the delta variances are
%
% $$ \Delta P_a = S P_a S' - \hat{P}_a, \quad  \Delta P_v = S P_v S' -
% \hat{P}_v, \quad \Delta P_w = S P_w S' $$
%
% There are several ways we can plot the sandpiles. When all the delta
% variances are positive, the sandpile should show the formal variance, and
% the deltas due to each component.  When all the deltas are negative, the
% sandpile should show the true variance, and negatives of the deltas.
% Otherwise, plot the components of the true variance as a positive
% sandpile, and the components of the formal variance as a negative
% sandpile, and relabel the negative y-axes to indicate this.

if nargout == 0,
    demomode = true;
else
    demomode = false;
end
if nargin == 0,
    % Demo/self-test mode: make up some data
    t = 0:100;
    lent = length(t);
    % Case 1: all delta variances positive:
    dVa = ones(1,1,lent);
    dVv = ones(1,1,lent);
    dVw = ones(1,1,lent);
    Va = shiftdim(exp(-t),-1);
    Vv = ones(1,1,lent) + 1;
    Vw = ones(1,1,lent) + 1;
    Vhata = Va - dVa;
    Vhatv = Vv - dVv;
    Vhatw = Vw - dVw;
    V = Va + Vv + Vw;
    Vhat = Vhata + Vhatv + Vhatw;
    figure % Open new figure window
    a = varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat);

    if demomode,
        disp('Hit any key to continue')
        pause
    else
        for k = length(a):-1:1,
            h1(k) = get(a(k));
        end
    end
    % Case 2: all delta variances negative:
    dVa = -ones(1,1,lent);
    dVv = -ones(1,1,lent);
    dVw = -ones(1,1,lent);
    Va = shiftdim(exp(-t),-1);
    Vv = ones(1,1,lent) + 1;
    Vw = ones(1,1,lent) + 1;
    Vhata = Va - dVa;
    Vhatv = Vv - dVv;
    Vhatw = Vw - dVw;
    V = Va + Vv + Vw;
    Vhat = Vhata + Vhatv + Vhatw;
    a = varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat);
    if demomode,
        disp('Hit any key to continue')
        pause
    else
        for k = length(a):-1:1,
            h2(k) = get(a(k));
        end
    end
    % Case 3: mixed-sign delta variances:
    dVa = -ones(1,1,lent);
    dVv = ones(1,1,lent);
    dVw = zeros(1,1,lent);
    Va = shiftdim(exp(-t),-1);
    Vv = ones(1,1,lent) + 1;
    Vw = ones(1,1,lent) + 1;
    Vhata = Va - dVa;
    Vhatv = Vv - dVv;
    Vhatw = Vw - dVw;
    V = Va + Vv + Vw;
    Vhat = Vhata + Vhatv + Vhatw;
    a = varpiles(t,dVa,dVv,dVw,Va,Vv,Vw,Vhata,Vhatv,Vhatw,V,Vhat);
    if demomode,
        disp('Hit any key to continue')
    else
        for k = length(a):-1:1,
            h3(k) = get(a(k));
        end
        close % Close current figure
    end
    if ~demomode,
        load varpiles_data
        try
            for k = length(h1):-1:1,
                fail1(k) = any(h1(k).YData - h1_save(k).YData);
            end
            for k = length(h2):-1:1,
                fail2(k) = any(h2(k).YData - h2_save(k).YData);
            end
            for k = length(h3):-1:1,
                fail3(k) = any(h3(k).YData - h3_save(k).YData);
            end
            fail = any([fail1,fail2,fail3]);
        catch %#ok<CTCH>
            fail = true;
        end
        varargout{1} = fail;
    end
    return
else
    t = varargin{1};
    t = t(:);
    dVa = squeeze(varargin{2});
    dVv = squeeze(varargin{3});
    dVw = squeeze(varargin{4});
    Va = squeeze(varargin{5});
    Vv = squeeze(varargin{6});
    Vw = squeeze(varargin{7});
    Vhata = squeeze(varargin{8});
    Vhatv = squeeze(varargin{9});
    Vhatw = squeeze(varargin{10});
    V = squeeze(varargin{11});
    Vhat = squeeze(varargin{12});
    
    if nargin > 12
        dVm = squeeze(varargin{13});
        Vm = squeeze(varargin{14});
        Vhatm = squeeze(varargin{15});
        deltastring = '\Delta Var Due to Man Exec';
        fullstring = 'Var Due to Man Exec';
    else
        dVm = [];
        Vm = [];
        Vhatm = [];
        deltastring = [];
        fullstring = [];
    end
end
tol = 1e-4*eps; % tolerance for plot type test
if all(all(dVa>tol)) && all(all(dVv>tol)) && all(all(dVw>tol)) && all(all(dVm>tol)),
    a = area(t, [Vhat, dVa, dVv, dVw, dVm]);
    legobjs = a;
    legstrs = { ...
        'Formal Variance', '\Delta Var Due to {\ita Priori}',...
        '\Delta Var Due to Meas Noise','\Delta Var Due to Proc Noise',...
        deltastring};
    if isempty(legstrs{5})
        legstrs(5) = [];
    end
elseif all(all(dVa<-tol)) && all(all(dVv<-tol)) && all(all(dVw<-tol)) && all(all(dVm<-tol)),
    a = area(t, [V, -dVa, -dVv, -dVw, -dVm]);
    legobjs = a;
    legstrs = { ...
        'True Variance', '\Delta Var Due to {\ita Priori}',...
        '\Delta Var Due to Meas Noise','\Delta Var Due to Proc Noise',...
        deltastring};
    if isempty(legstrs{5})
        legstrs(5) = [];
    end
else
    a = area(t,[-Vhata,-Vhatv,-Vhatw,-Vhatm]);
    set(a(1),'facecolor',[0 .75 1])
    set(a(2),'facecolor',[1 .75 0])
    set(a(3),'facecolor',[.5 .05 .03])
    if ~isempty(Vhatm)
        set(a(4),'facecolor',[.5 1 .2])
    end
    hold on
    a0 = a;
    a = area(t,[Va,Vv,Vw,Vm]);
    set(a(1),'facecolor',[0 .75 1])
    set(a(2),'facecolor',[1 .75 0])
    set(a(3),'facecolor',[.5 .05 .03])
    if ~isempty(Vhatm)
        set(a(4),'facecolor',[.5 1 .2])
    end
    hold off
    fixYLabel;
    axis manual
    set(get(gca,'ylabel'),'Units','Normalized')
    ylpos = get(get(gca,'ylabel'),'position');
    text(ylpos(1),0.25,'Formal Variance','HorizontalAl','Center',...
        'VerticalAl','Bot','Rotation',90,'Units','Normalized')
    text(ylpos(1),0.75,'True Variance','HorizontalAl','Center',...
        'VerticalAl','Bot','Rotation',90,'Units','Normalized')
    legobjs = a;
    legstrs = { ...
        'Var Due to {\ita Priori}','Var Due to Meas Noise',...
        'Var Due to Proc Noise',fullstring};
    a = [a0,a];
    if isempty(legstrs{4})
        legstrs(4) = [];
    end
    z_handle = zoom(gcf);
    set(z_handle,'ActionPostCallback',@fixYLabel)
end
legend(legobjs,legstrs)
title('Variance Sandpile')
zoom yon
varargout{1} = a;

end

function fixYLabel(~,~)

set(gca,'YticklabelMode','auto')
lim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-lim lim])
set(gca,'yticklabel',num2str(abs(get(gca,'ytick')')))

end