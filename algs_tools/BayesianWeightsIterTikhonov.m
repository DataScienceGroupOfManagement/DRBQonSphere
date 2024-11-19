function [alpha_seq] = BayesianWeightsIterTikhonov(xtr,  ytr, KerPara,  lambda_seq, iterMax )
    Ntr= size(xtr,1);
    alpha_seq = zeros(Ntr, length(lambda_seq));
    num_lambdas = length(lambda_seq);
    Ktr = KernelComputation(xtr, xtr, KerPara); 
    Ktr = (Ktr + Ktr')/2;
    if num_lambdas>1
        Ktr = (Ktr+Ktr')/2;
        [U, S] = eig(Ktr); 
        diagS = diag(S);
    end
    for k = 1:1: num_lambdas
        lambda_value = lambda_seq(k);
        if num_lambdas == 1
            Kinv = inv(Ktr+Ntr*lambda_value*eye(Ntr));
        else
            Kinv = U*diag(1./(diagS+Ntr*lambda_value))*U';
        end
        alpha = zeros(Ntr,1);
        for iter = 1:iterMax
            alpha_new = Kinv*(ytr+Ntr*lambda_value*alpha);
            alpha = alpha_new;
        end
        alpha_seq(:, k) = alpha;
    end
end