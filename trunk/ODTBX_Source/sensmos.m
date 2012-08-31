function varargout = sensmos(t,Sig_sa)
% Sensitivity Mosaic Plot: A Checkerboard plot of the sensitivity matrix.

% SENSMOS(T,SIG_SA) plots a "sensitivity mosaic," which is a checkerboard
% plot of the sensitivity matrix, SIG_SA.  The sensitivity is plotted at
% the final time in the vector T.  A slider allows the user to scroll to
% other times.  Warm colors indicate a strong sensitivity, while cool
% colors indicate a weak sensitivity.
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
%
%  REVISION HISTORY
%   Author      		    Date         	Comment
%   Unknown                Unknown          Created
%   Ravi Mathur          08/27/2012         Extracted regression test to
%                                           its own file (see regression
%                                           testing framework)

Sig_sa = 10*log10(abs(eps+Sig_sa));
[ns,n,el] = size(Sig_sa);
pcolor(eye(ns+1,ns)*Sig_sa(:,:,end)*eye(n,n+1))
th = text(0.5,-.1,datestr(t(end),'HH:MM'),...
    'units','normalized','horizontalAlign','center',...
    'tag','mostime');
set(gca,'xtick',1:n,'ytick',1:ns)
xlabel('{\it A Priori} State Index')
sa = uicontrol(gcf,'style','slider',...
    'max',el,'min',1,...
    'value',el,...
    'sliderstep',[1 10]./(el-1),...
    'units','normalized','position',...
    get(gca,'position')*[1 0 0 0;0 1 0 -.1;0 0 1 0;0 0 0 .1]');
    function ca(~,~)
        j = get(sa,'value');
        pcolor(eye(ns+1,ns)*Sig_sa(:,:,j)*eye(n,n+1));
        set(th,'string',datestr(t(j),'HH:MM'))
    end
set(sa,'callback', @ca)
ylabel('Solve-For State Index at Slider Time')
set(gca,'ydir','rev','xaxisloc','top')
axis tight
hc = colorbar;
set(get(hc,'title'),'string','dB')
title('Logarithmic Sensitivity Mosaic')
hold on
% This will scale colorbar to cover the whole range
for i = 1:el
    pcolor(eye(ns+1,ns)*Sig_sa(:,:,i)*eye(n,n+1))
    set(sa,'value',i)
end
end