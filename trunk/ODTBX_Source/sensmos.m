function sensmos(t,Sig_sa)
% SENSMOS(T,SIG_SA) plots a "sensitivity mosaic," which is a checkerboard
% plot of the sensitivity matrix, SIG_SA.  The sensitivity is plotted at
% the final time in the vector T.  A slider allows the user to scroll to
% other times.  Warm colors indicate a strong sensitivity, while cool
% colors indicate a weak sensitivity.

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
        i = get(sa,'value');
        pcolor(eye(ns+1,ns)*Sig_sa(:,:,i)*eye(n,n+1));
        set(th,'string',datestr(t(i),'HH:MM'))
    end
set(sa,'callback', @ca)
ylabel('Solve-For State Index at Slider Time')
set(gca,'ydir','rev','xaxisloc','top')
axis tight
hc = colorbar;
set(get(hc,'title'),'string','dB')
title('Logarithmic Sensitivity Mosaic')
hold on
end