function [alpha_seq ] = BayesianWeightsIterLandweber(xtr,  ytr, KerPara,  tau_seq, iterMax)
    Ntr= size(xtr,1);
    alpha_seq = zeros(Ntr, length(tau_seq));
    num_taus = length(tau_seq);
    Ktr = KernelComputation(xtr, xtr, KerPara); % the kernel matrix
    Ktr = (Ktr + Ktr')/2;
    for k = 1:1: num_taus
        tau = tau_seq(k);
        alpha = zeros(Ntr,1);
        for iter = 1:1:iterMax
            alpha_new = alpha+(ytr-Ktr*alpha)*tau/Ntr;
            alpha = alpha_new;
        end
        alpha_seq(:, k) = alpha;
    end
end