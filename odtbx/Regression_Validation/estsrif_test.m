function failed = estsrif_test
%
% estsrif_test Regression testing and demo for estsrif.m
% See also estsrif.m
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)
%
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
%   Ravi Mathur             08/28/2012      Extract regression test from
%                                           original estsrif.m
%
%
%
% Note that the brunt of the regression testing still occurs in
% estsrif.m. It is deeply embedded and will require considerable effort
% to extract all of it into estsrif_test.m.

totaltests = 3;
disp('estsrif_test: regression testing estsrif...')

% Run all the test cases
for k = 1:totaltests,
    disp(['Case ',num2str(k),'...'])
    fail(k,:) = estsrif(k); %#ok<AGROW>
end

% If system supports parallel processing run all
% the test cases again in parallel
testparallel = 0;
try % See if system supports parallel operation
    if matlabpool('size') == 0,
        matlabpool('open');
        testparallel = 1;
    else
        disp('Skipping parallel test.');
        disp('Self test already running in parallel environment.');
    end
catch %#ok<CTCH>
    disp('Skipping parallel test.');
    disp('No support for parallel processing.');
end
if testparallel,
    for k = 1:totaltests,
        disp(['Parallel Case ',num2str(k),'...'])
        fail(k+totaltests,:) = estsrif(k); %#ok<AGROW>
    end
    matlabpool('close');
end

failed = any(any(fail));