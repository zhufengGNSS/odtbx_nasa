function fixzoom
% This fixes the zoom function so that zero stays in view when vertically
% zooming.
h = zoom(gcf);
set(h,'Motion','vertical','Enable','on','ActionPostCallback',@keepzero,...
    'ActionPreCallback',@saveylims)

function saveylims(obj,evd) %#ok<INUSL>
% This is needed so postcallback can know what to do
set(evd.Axes,'UserData',get(evd.Axes,'Ylim'))

function keepzero(obj,evd) %#ok<INUSL>
% This callback puts zero back in view after zooming
oldLim = get(evd.Axes,'UserData');
newLim = get(evd.Axes,'YLim');
if oldLim(1) == 0,
    set(evd.Axes,'Ylim',[0, newLim(2)])
else
    maxLim = max(abs(newLim));
    set(evd.Axes,'Ylim',[-maxLim, maxLim])
end