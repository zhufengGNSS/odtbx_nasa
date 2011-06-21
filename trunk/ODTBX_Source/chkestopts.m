function estopts = chkestopts(options,ncases,m)
% CHKESTOPTS Check some of the estimator options for consistency
%
% Checks the consistency of the specifications for the Monte Carlo seed and
% measurement editing options
% 
%   estopts = chkestopts(options,ncases,m)
%
% INPUT: 
%  options   Options structure from odtbxOptions
%  ncases    Number of Monte Carlo cases
%  m         Number of measurements (If not specified, measurement editing
%            options are not checked.)
% 
% OUTPUT:
%  estops   Checked options
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

% Created by S. Hur-Diaz
% Emergent Space Technologies
% 7/30/2010

%% Check Monte Carlo seed
monteseed     = getOdtbxOptions(options, 'MonteCarloSeed', NaN);
monteseed_use = NaN(1,ncases); % Pre-allocate array

if(~isnan(monteseed))
    if((length(monteseed) ~= ncases) && (length(monteseed) ~= 1))
        error('Number of Monte Carlo Seeds specified does not match the number of cases.');
    elseif((ncases > 1) && (length(monteseed) == 1))
        for bb=1:ncases
            monteseed_use(bb) = monteseed + (bb-1);
        end
    else
        monteseed_use = monteseed;
    end
else
    for bb=1:ncases
        monteseed_use(bb) = monteseed;
    end
end
estopts.monteseed = monteseed_use;

%% Check measurement editing options
if nargin > 2
    
    editvec = getOdtbxOptions(options, 'EditVector', 0);
    eratio = getOdtbxOptions(options, 'EditRatio', []);
    
    if isempty(eratio)
        
        warning('CHKESTOPTS:EditRatio',' EditRatio will be set to the default value of 9');
        eratio = 9*ones((m-1)*(1-editvec)+1,1);
        
    elseif editvec == 1 && length(eratio) ~= 1
        
        error('dimERATIO: dim(ERATIO)~= 1 but EDITVEC is set to 1.');
        
    elseif length(eratio) < m && editvec == 0
        
        warning('CHKESTOPTS:dimERATIO',...
            ['dim(ERATIO)~=dim(Y); ERATIO(1) used for all.\n To suppress ', ...
            'this message, use warning(''off'', ''KALMUP:dimERATIO'')']);
        eratio = repmat(eratio(1),m,1);
        
    end
    
    eflag  = getOdtbxOptions(options, 'EditFlag', []);
    
    if isempty(eflag)
        
        warning('CHKESTOPTS:EditFlag', 'EditFlag will be set to the default value of 1');
        eflag = ones((m-1)*(1-editvec)+1,1);
        
    elseif editvec == 1 && length(eflag) ~= 1
        
        error('dimEFLAG: dim(EFLAG)~= 1 but EDITVEC is set to 1.');
        
    elseif length(eflag) < m && editvec == 0
        
        warning('CHKESTOPTS:dimEFLAG',...
            ['dim(EFLAG)~=dim(Y); EFLAG(1) used for all.\n To suppress ', ...
            'this message, use warning(''off'', ''KALMUP:dimELAG'')']);
        eflag = repmat(eflag(1),m,1);
        
    end
    
    estopts.eratio = eratio;
    estopts.eflag = eflag;
    estopts.editvec = editvec;
end
