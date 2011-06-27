function closegmat()

% CLOSEGMAT GMAT engine shutdown routine
%  closegmat Closes the GMAT interface and unloads the core library from 
%  the workspace. 
%
%   keyword: GMAT Integrators Forces 
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

% Unload the library
fprintf(1,'Unloading libCInterface...\n');
unloadlibrary('libCInterface');

if libisloaded('libCInterface')
   error('Unloading failed\n');
else
   fprintf(1,'Unloaded\n');
end    

end %function
