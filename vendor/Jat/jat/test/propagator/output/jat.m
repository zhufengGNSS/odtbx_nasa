function X = jat(file)
% jat
% simple function to read in the JAT output into the Matlab workspace

X = dlmread(file,'\t',0,0);

% max2 = max(max(X(1:num,2)),max(X(1:num,3)));
% min2 = min(min(X(1:num,2)),min(X(1:num,3)));
% max3 = max(max2,max(X(1:num,4)))
% min3 = min(min2,min(X(1:num,4)))
%figure(1);
%plot3(X(1:num,2),X(1:num,3),X(1:num,4));
% xlim([min3 max3]);
% ylim([min3 max3]);
% zlim([min3 max3]);
%xlabel('x');
%ylabel('y');
%zlabel('z');
%figure(2);
%plot(X(1:num,2),X(1:num,3));
% xlim([min2 max2]);
% ylim([min2 max2]);
%xlabel('x');
%ylabel('y');