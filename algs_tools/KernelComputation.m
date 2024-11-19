% calculate the kernel matrix.
% X is a mxd matrix，Y is a nxd matrix，each row is a data sample of d-dim
% type = 1 : K(x,y) = 1+min(x,y);
%      = 2 : K(x,y) = h(||x-y||_2),
%                       /  (1-t)^4*(4*t+1),    if  0<t<=1
%            其中h(t) = |
%                       \   0                  if  t>1
%      = 4 : K(x,y) = exp{-||x-y||_2^2/(2*Sigma^2)}

function K = KernelComputation(X, Y, KerPara)
m = size(X,1);
n = size(Y,1);
type = KerPara.KernelType;
if isfield(KerPara,'para')
    para = KerPara.para;
end
switch type
    case 1
        XX = X*ones(1,n);
        YY = ones(m,1)*Y';
        K = ones(m,n) + (XX + YY)/2 - abs(XX - YY)/2;
    case 2 % the integral value is  pi/7
        KK = pdist2(X,Y,'euclidean');
        temp_index = KK>=1;
        K = (1 - KK).^4.*(4*KK + 1); % function of C^{2}
        K(temp_index) = 0;
    case 3
        KK = pdist2(X,Y,'euclidean');
        temp_index = KK>=1;
        K = (1 - KK).^6.*(35*KK.^2+18*KK+3); % function of C^{4}
        K(temp_index) = 0;
    case 4    % gaussian kernel
        Sigma = para;
        KK = pdist2(X,Y,'euclidean');
        KK = KK/Sigma;
        K = exp((-KK.^2)/2);
    case 5
        k = para(1);
        delta = para(2);
        r = pdist2(X, Y, 'euclidean');  % n.points X n.center
        K = Wendland_r(r, k, delta);
    case 6  %  MQ (multiquadric) Kernel
        Sigma = para;
        EpsilonSquare = 1/(2*Sigma*Sigma);
        KK = pdist2(X,Y,'euclidean');
        K = sqrt(1+EpsilonSquare.*KK.^2);
    case 7  % IMQ (inverse multiquadric)
        Sigma = para;
        EpsilonSquare = 1/(2*Sigma*Sigma);
        KK = pdist2(X,Y,'euclidean');
        K = 1 ./sqrt(1+EpsilonSquare.*KK.^2);
    case 8 % (1-r)^8 * (32r^3+25r^2+8r+1)  function of C^{6}
        KK = pdist2(X,Y,'euclidean');
        temp_index = KK>=1;
        K = (1 - KK).^8.*(32*KK.^3 + 25*KK.^2 +8*KK+1);
        K(temp_index) = 0;
    case 9 % (1-2r)^4  * (8r+1) , 0<=r<=1
        KK = pdist2(X,Y,'euclidean');
        temp_index = KK>=1;
        K = (1 - 2*KK).^4.*(4*KK*2 + 1);
        K(temp_index) = 0;
end
