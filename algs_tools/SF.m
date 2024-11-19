function x = SF(L)
% [w,y] = SD(L,d)
% load the symmetric spherical design of Rob

% Inputs:
% L -- degree for which the SD quadrature is exact
% Outputs:
% y -- nodes of quadrature rule, size(y)=[Num.points, 3]

loadfpath = 'Points/SFD/';

if L<10
    Ltxt = ['00' num2str(L)];
elseif L<100
    Ltxt = ['0' num2str(L)];
else
    Ltxt = num2str(L);
end
ld = [loadfpath 'sf' Ltxt '.mat'];
load(ld,'data');
x = data;