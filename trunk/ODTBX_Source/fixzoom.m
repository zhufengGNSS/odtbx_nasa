function fixzoom
% Fixes the zoom function so that zero stays in view when vertically zooming.

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