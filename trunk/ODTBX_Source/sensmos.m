function varargout = sensmos(t,Sig_sa)
% SENSMOS(T,SIG_SA) plots a "sensitivity mosaic," which is a checkerboard
% plot of the sensitivity matrix, SIG_SA.  The sensitivity is plotted at
% the final time in the vector T.  A slider allows the user to scroll to
% other times.  Warm colors indicate a strong sensitivity, while cool
% colors indicate a weak sensitivity.

if nargin == 0 % demo/test mode
    t = (1:12)./24;
    for i = length(t):-1:1
        S(:,:,i) = t(i)*ones(6,8);
    end
    try
        sensmos(t,S)
        fail = false;
    catch %#ok<CTCH>
        fail = true;
    end
    varargout{1} = fail;
    return
end
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