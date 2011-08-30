function [retval] = opengmat;

% OPENGMAT GMAT engine startup routine
%  [retval] = opengmat() Loads the GMAT interface and core library 
%
%   OUTPUTS
%   retval          (1x1)       Integer return code indicating the status
%                               of the GMAT startup process.
%
%   keyword: GMAT 
%
%
% (This file is part of GMAT, The General Mission Analysis Tool, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

% GMAT: General Mission Analysis Tool
%
% Copyright (c) 2002-2011 United States Government as represented by the
% Administrator of The National Aeronautics and Space Administration.
% All Other Rights Reserved.
%
% Developed jointly by NASA/GSFC and Thinking Systems, Inc. under NASA
% Prime Contract NNG10CP02C, Task Order 28.
%
% Author: Darrel J. Conway, Thinking Systems, Inc.
% Created: 2011/05/17

% Load the library
if ~(libisloaded('libCInterface'))
    [notfound, warnings] = loadlibrary('libCInterface', ...
        @interfacewrapper);
    
    if (~isempty(notfound)) || (~libisloaded('libCInterface'))
       error('Unable to load the libCInterface shared library.  Warnings were: %s\n', warnings);
    end
end

% Start the GMAT engine
retval = calllib('libCInterface','StartGmat');
if retval ~= 0
    error('GMAT failed to start.  This may be caused by missing or misconfigured files; check the gmat_startup_file.\n');
end
msg = calllib('libCInterface','getLastMessage');

disp(msg);

end % function
