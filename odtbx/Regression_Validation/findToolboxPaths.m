function [toolboxPaths, otherPaths] = findToolboxPaths()
% Parses the MATLAB search paths and tries to determine which paths
% are for optional MATLAB toolboxes.  The toolboxPaths is a cell array
% of strings that are MATLAB toolbox paths.  Some paths with "toolbox"
% in them are ignored.  All paths that aren't directly MATLAB-related
% are placed in otherPaths, a cell array of strings.  The standard
% MATLAB paths under the toolbox directory are ignored.  Case-insensitive
% path comparisons are used.
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

toolboxPaths = {};
otherPaths = {};

fs = filesep;

% all MATLAB toolboxes are in the MATLAB directory, with a release number R20*
% they also live in the toolbox directory
ML{1} = lower(strcat(fs,'MATLAB',fs,'R20')); 
ML{2} = lower(strcat(fs,'toolobx',fs)); 

% Paths under MATLAB to keep:
keep = {lower(strcat(fs,'toolbox',fs,'matlab',fs)) ...
           lower(strcat(fs,'toolbox',fs,'local'))...
           lower(strcat(fs,'toolbox',fs,'shared'))...
           lower(strcat(fs,'work')) };

% the paths split into a single cell of n strings, use a large buffer
% for boxes with lots of toolboxes
p = textscan(path,'%s','Delimiter',pathsep,'BufSize',20480);

for i = 1:length(p{1})
    pstring = char(p{1}(i));
    if isempty(findstr(lower(pstring),ML{1})) && isempty(findstr(lower(pstring),ML{2}))
        otherPaths{end+1} = pstring;
    else
        keepit = 0;
        for j = 1:length(keep)
            if ~isempty(findstr(lower(pstring),char(keep{j})))
                keepit = 1;
                break;
            end
        end
        if keepit == 0
            toolboxPaths{end+1} = pstring;
        end
    end
end
