function [xdot A Q] = dualIADyn(t,x,~)

lent = length(t);

I = eye(6,6);
O = zeros(6,6);

xc = x(1:6,:);
xt = x(7:12,:);

[xdotC Ac Qc] = r2bp(t,xc);
[xdotT At Qt] = r2bp(t,xt);

xdot = [xdotC;
        xdotT];

A = nan(12,12,lent);
Q = nan(12,12,lent);
for i=1:lent
    A(:,:,i) = [Ac O;
                O At];
    Q(:,:,i) = blkdiag(Qc,Qt);
end
    
end