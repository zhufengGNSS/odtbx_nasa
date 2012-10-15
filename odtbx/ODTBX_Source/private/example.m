function y=Example(x)

%EXAMPLE This is an example function header. (must be on first line only for contents.m.)
%   EXAMPLE solves the ... problem (Descriptive text of how the function works.)
%  
%   INPUTS
%   VARIABLE    SIZE    DESCRIPTION (Optional/Default)
%      X        (3xN)     Input to function
%
%   OUTPUTS 
%      Y        (3XN)     Output from function
%                                       
%   *Show different variations of the function: bintprog example from Optimization Toolbox 
%
%   DEFAULT
%   Y = EXAMPLE(X)  Add this since it would allow us to see the
%   arguments being passed in right away.
%
%   OPTIONAL
%   X = BINTPROG(f,A,b) solves the problem min f'*X subject to the linear 
%   inequalities A*X <= b, where the elements of X are binary integers.
%
%   *Give an example of how the function works:
%
%   EXAMPLE
%     f = [-9; -5; -6; -4]; 
%     A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
%     b = [9; 1; 0; 0];
%     X = bintprog(f,A,b) 
%
%    keyword: Time, GPS, Date, Coordinate Transformations, Plotting Utilities/Graphics, File I/O, 
%
%   * Be careful not to add any uncommented spaces in the comment block. If you do,
%   any comments below the initial comment block will be ignored by the
%   command window and will not be displayed when you type help.
%   * This header can be used for any file type just edit it to be
%   applicable for your file type. 
%   * See also must be last line or else anything in caps will be
%   hyperlinked by the See also command. The link will not work. 
%   See also BINTPROG
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

%   REVISION HISTORY
%   Author      Date         Comment
%               (MM/DD/YYYY)
%   Bo Naasz     09/01/2005   Original
