function y = sinc_multicentre(xyz, xc, c, w)
% Y = sinc_multicentre(xyz,k,delta,xc,w)
% is a linear combination of sinc functions on sphere with centres xc and
% weigts w.
%
% Inputs:
% xyz -- points to evaluate, size(xyz,1)=Num. Points, size(xyz,2)=3
% k -- type of Wendland; smoothness = k + 3/2
% delta -- Scaling factor delta > 0 (Default delta = 1)
% xc -- set of centres of the Wendland functions;
% size(xc) = [Num. centres, Dim. sphere +1]
% w -- set of weights; row vector

d = size(xyz, 2);

if nargin < 2
    xc = [eye(d) -eye(d)]';
end

if nargin < 3
    c = [1 1];
end

if nargin < 4
    w = ones(size(xc, 1),1);
end

r = c(1) * acos(xyz * xc'); % n.points X n.center
Y = sin(c(2)*pi*r)./(r);

y = Y * w;
