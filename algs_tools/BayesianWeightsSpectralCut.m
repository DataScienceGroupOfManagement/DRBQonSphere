function [alpha_seq] = BayesianWeightsSpectralCut( xtr,  ytr, KerPara,  lambda_seq )
    Ntr= size(xtr,1);
    alpha_seq = zeros(Ntr, length(lambda_seq));
    Ktr = KernelComputation(xtr, xtr, KerPara);
    Ktr = (Ktr + Ktr')/2;
    [U,S,V] = svd(Ktr);
    U= real(U);
    S= real(S);
    tempS = 1./diag(S);
    for k = 1:1:length(lambda_seq)
        lambda_value = lambda_seq(k);
        temp = 1./diag(S);
        temp(tempS>(1/lambda_value)) = 0;
        diagSinv = diag(temp);
        KlambdaInv = U * diagSinv * V';
        alpha = KlambdaInv*ytr;
        alpha_seq(:, k) = alpha;
    end
end