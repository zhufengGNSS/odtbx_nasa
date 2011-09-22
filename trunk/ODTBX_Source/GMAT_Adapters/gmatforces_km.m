function [xDot,A,Q] = gmatforces_km(t,x,options)

% GMATFORCES_KM  Returns derivatives from GMAT ODE models (state in km, sec).
%
%   xDot = gmatforces_km(t,x,~) returns the derivatives of the GMAT
%   force models specified by the given modelid.  If modelid is omitted,
%   the first odemodel found will be used.  These models are created in a
%   GMAT script that must be loaded using preparescript.
%
%   [xDot,A] = gmatforces_km(t,x,~) also returns the Jacobian
%   (km, sec)
%
%   [xd,A,Q] = gmatforces_km(t,x,~) also returns a process noise
%   power spectral density based on a random walk process noise input.
%   At present the magnitude of the PSD is hardcoded to be 1e-9 (km, sec).
%
%   INPUTS
%   VARIABLE        SIZE    	DESCRIPTION (Optional/Default)
%       t           (1xn)       Time since start of simulation (secs)
%       x           (6xn)       Input state (x,y,z position in km
%                               followed by x,y,z velocity in km/s)
%    modelid        (1x1)       Integer id for forcemodel (optional)
%
%   OUTPUTS
%       xDot        (6xn)       Derivatives of state (km, seconds)
%       A           (6x6xn)     State transition matrix (km, seconds)
%       Q           (6x6xn)     process noise power spectral density
%
%   keyword: GMAT Forces
%   See also startgmat, closegmat, preparescript

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
% Modified: 2011/09/06 Russell Carpenter, NASA GSFC

if ~(libisloaded('libCInterface'))
    error('Please load GMAT using the command startgmat before calling gmatforces');
end

[nx,nt]=size(x);
if nx ~= 6
    error('gmatforces_km expects a state size of 6, given %d', nx);
end

xDot=zeros(nx,nt);
A = zeros(nx,nx,nt);

if (nt > 1)
    fprintf(1,'Warning: The interface is not yet tested for more than one spacecraft!!!');
end

t0 = getOdtbxOptions(options,'epoch',0)/86400;
modelid = getOdtbxOptions(options,'modelid',[]);

if ~isempty(modelid)
    setodemodel(modelid);
end

for i=1:nt
    % Set the state.  Note that GMAT needs an a.1 epoch.
    epoch = t0 + 21545.0 + t(i) / 86400.0;
    statep = libpointer('doublePtr', x(:,i));
    retval = calllib('libCInterface','SetState',epoch,statep,nx);
    if (retval < 0)
        fprintf(1,'Failure in libCInterface SetState()\n');
        errmsg = calllib('libCInterface','getLastMessage');
        fprintf(1,'The error message was: %s\n',errmsg);
        return;
    end
    
    % Call GetDerivatives
    dim = 0;  % Note, this value won't change
    dimptr = libpointer('int32Ptr', dim);
    derivsptr = calllib('libCInterface','GetDerivatives', 0.0, 1, dimptr);
    dimval = get(dimptr,'value');
    if dimval == 0
        fprintf(1,'Failure in libCInterface GetDerivatives() dimval:\n');
        errmsg = calllib('libCInterface','getLastMessage');
        fprintf(1,'The error message was: %s\n',errmsg);
        return;
    end
    
    % Set the return data
    setdatatype(derivsptr, 'doublePtr', dimval, 1);
    derivatives = get(derivsptr,'Value');
    xDot(:,i) = derivatives(1:6);
    
    if nargout > 1
        for m=1:6
            for n=1:6
                A(m,n,i) = derivatives(m*6+n);
            end
        end
    end
end

if nargout == 3,
    Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9].^2),[1 1 nt]);
end

end % function
