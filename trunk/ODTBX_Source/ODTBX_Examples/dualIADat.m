function [y H R] = dualIADat(t,x,options)

sig = options.sig;

lent = length(t);
I = eye(3,3);
O = zeros(3,3);

xc = x(1:3,:);
xt = x(7:9,:);

q = repmat([0 0 0 1]',1,lent);
[yg Hg Rg] = gpsmeas(t,x(1:6,:),options,q);

n = size(yg,1);

yp = xt-xc;

y = [yp;yg];

for i=lent:-1:1
    H(:,:,i) = [-I O I O; Hg(:,:,i) zeros(n,6)];
    R(:,:,i) = blkdiag(diag(sig.^2),Rg(:,:,i));
end

end