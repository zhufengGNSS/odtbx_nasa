function [indA, matchA] = compValsTol(A, B, tol)
% Compares values using a tolerance.
%
% Compares the values of vector A to those in vector B using a tolerance,
% tol.  The indices of vector A that match one or more values in vector B
% are returned in indA.  The tolerance for a match is applied as:
%
% B(i)-tol <= A(j) <= B(i)+tol
%
% Both A and B are assumed to be 1xN row vectors of unique doubles.  If tol
% is not supplied it is assumed to be 0.0.
%
% Note, this function uses vectorized math and has several orders of
% magnitude performance improvement over brute-force loops or some other 
% diff/sort approaches.
%
% Usage:
%   [indA, matchA] = compValsTol(A, B, tol);
%
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   A               double  1xN     Array of values to compare with B
%   B               double  1xM     Array of values being searched and
%                                   compared
%
%   INPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION (Optional/Default)
%   indA            double  1xX     The indices of A that compare within 
%                                   tolerance to one or more values in B, 
%                                   could be empty.
%   matchA          double  1xX     The values of A that compare within 
%                                   tolerance to one or more values in B, 
%                                   could be empty.
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

indA = [];

if isempty(A) && isempty(B)
    return;
end

if ~exist('tol','var') || isempty(tol)
    tol = 0;
end

ab = [A B]; % combine the data

% combined array to hold the data and a source tag
abc = zeros(2,length(ab)); % (1,:)=data, (2,:)=source tag

% sort combined, get indices
[abc(1,:),absind] = sort(ab);

% label the source
sortedtinds = absind > length(A);
abc(2,:) = 1; % 1=A
abc(2,sortedtinds) = 2; % 2=B

% differences in source
ds = diff(abc(2,:),1,2); % numerical

% differences in value, comapred to the tol
ltol = abs(diff(abc(1,:),1,2)) <= tol; % logical

% Note, this logic works on the diff arrays, so the indexing is off by 1
% compared to abc indices.
pass1 = (ds == 1) & ltol; % check B with A to its left
indabc1 = find(pass1); % these A locations use the diff index (one shorter than abs index)

pass2 = (ds == -1) & ltol; % check B with A to its right
indabc2 = find(pass2)+1; % these A locations are one to the right of the diff index

indabc = sort([indabc1 indabc2]);

indA = absind(indabc);
if nargout > 1
    matchA = A(indA);
end
