function X = unifsphere(n, d)
    % n - number of samples
    % d- dimensinality
    sigma = eye(d);
    mu = zeros(n,d);
    M = mvnrnd(mu,sigma);
    M_norm = sqrt(sum(M.^2,2));
    X = M./repmat(M_norm,1,d);
end  