function failed = convertTime_test()
% Regression Test Case
% Function(s) convertTime
%
% This function tests the time conversion function to/from GPS, UTC, TT,
% TDB, TAI
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

% Modified by 
%     Keith Speckman      3/18/2008      Modified to reflect code change of GPS origin date to Jun 6, 1980
%                                        Corrected "any(any(..." line which always resutled in a failure
%     Allen Brown         10/15/2009     Corrected test to properly apply
%                                        GPS epoch conversion to test inputs and
%                                        outputs.

failed  = 0;
tol     = 1e-9; % 0.0864 ms (in datenum scale)
epoch   = datenum('June 1, 2007'); % datenum epoch (not GPS epoch)

% Output values (without June 6, 1980 GPS correction) have been commented out.
%                     to GPS                      to UTC                    to TT                     to TDB                   to TAI
% targetOutput = [7.331940000000000e+005    7.331939998379629e+005    7.331940005924074e+005    7.331940005924179e+005    7.331940002199074e+005
%                 7.331940001620371e+005    7.331940000000000e+005    7.331940007544444e+005    7.331940007544550e+005    7.331940003819445e+005
%                 7.331939994075926e+005    7.331939992455556e+005    7.331940000000000e+005    7.331940000000105e+005    7.331939996275000e+005
%                 7.331939994075819e+005    7.331939992455449e+005    7.331939999999893e+005    7.331940000000000e+005    7.331939996274893e+005
%                 7.331939997800926e+005    7.331939996180555e+005    7.331940003725000e+005    7.331940003725105e+005    7.331940000000000e+005 ];

% Output values with epoch corrections included.
% Those in the "to GPS" times are in the GPS epoch.
% All other "to" times are in the datenum epoch.
%                     to GPS                      to UTC                    to TT                     to TDB                   to TAI
targetOutput = [1.000800000000000e+04     7.331939998379629e+05     7.331940005924074e+05     7.331940005924179e+05     7.331940002199074e+05    % from GPS
                1.000800016203709e+04     7.331940000000000e+05     7.331940007544444e+05     7.331940007544550e+05     7.331940003819445e+05    % from UTC
                1.000799940759258e+04     7.331939992455556e+05     7.331940000000000e+05     7.331940000000105e+05     7.331939996275000e+05    % from TT
                1.000799940758187e+04     7.331939992455449e+05     7.331939999999893e+05     7.331940000000000e+05     7.331939996274892e+05    % from TDB
                1.000799978009262e+04     7.331939996180555e+05     7.331940003725000e+05     7.331940003725105e+05     7.331940000000000e+05 ]; % from TAI

formats = {'GPS', 'UTC', 'TT', 'TDB', 'TAI'};
n       = length(formats); 

output  = zeros(n,n);

for i=1:n
    for j=1:n
        if i == 1
            % from GPS timescale, so convert input epoch from datenum to
            % GPS epoch:
            epochin = epoch - datenum('Jan 6 1980');
        else
            % not GPS input, so don't alter the input epoch
            epochin = epoch;
        end
        output(i,j) = convertTime( formats{j}, formats{i}, epochin );
    end
end

diff = output-targetOutput;

if( any( any( abs(diff) > tol ) ) )
    failed = 1;
end

    epoch1 = datenum('June 1, 2007');
    epoch2 = datenum('June 1, 2008');
    epoch3 = datenum('June 1, 2009');

    epoch_vec = [epoch1 epoch2 epoch3];
    
    size_check = size(epoch_vec);
    
    formats = {'GPS', 'UTC', 'TT', 'TDB', 'TAI'};
    
    format_pairs{1,:} = {'UTC', 'TT',  'TDB', 'TAI'};
    format_pairs{2,:} = {'GPS', 'TT',  'TDB', 'TAI'};
    format_pairs{3,:} = {'GPS', 'UTC', 'TDB', 'TAI'};
    format_pairs{4,:} = {'GPS', 'UTC', 'TT',  'TAI'};
    format_pairs{5,:} = {'GPS', 'UTC', 'TT',  'TDB'};

    m = length(formats);
    n = length(format_pairs{1});
    
    k = 0;
    
    ctr = 0;
    
    for a=1:m
       
        for b=1:n
            
            res = convertTime(formats{a}, format_pairs{a}{b}, epoch_vec);
           
            size_res = size(res);
            
            if(size_res == size_check)
                % Ok, no problems here.
            else
                % Oops, we have a problem here.
                k = k + 1;
               
                problems{k,:} = {formats{a}, format_pairs{a}{b}}; %#ok<AGROW>
                
                prob_sizes(k,:) = size_res; %#ok<AGROW>
                
            end
 
            ctr = ctr + 1;
            
        end
        
    end
       
    if(k == 0)
       
        % No problem.
        
    else
        
        failed = 1;
        
        fprintf(1, '\n');
        
        fprintf(1, '%d problem(s) found out of %d possible permutations.\n', k, ctr);
        
        fprintf(1, '\n');
        
        for a=1:k
            fprintf(1, '%d: %s to %s, result size = [%d %d], should be [%d %d].\n', a, problems{a}{1}, problems{a}{2}, prob_sizes(a,1), prob_sizes(a,2), size_check(1), size_check(2));
        end
        
        fprintf(1, '\n');
        
    end
    
    
