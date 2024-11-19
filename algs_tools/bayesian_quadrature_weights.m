function [weights_seq] = bayesian_quadrature_weights(K, Vphi0, lambda_seq)
% BAYESIAN_QUADRATURE_WEIGHTS  
% K is the kernel matrix . N (points) x N (points)
% Vphi0 = (phi0, ..., phi0)^T is the vector of  kernel mean embedding U(xi).  N (points) x 1
% lambda_seq is the list of Regularization Parameter of KRR. (Carefully choose this)

% weights_seq : the t-th column is the bayesian quadrature weights according to
% the t-th parameter of lambda_seq

if nargin < 3
    lambda_seq = 1; 
end

KlambdaInv_seq = K_add_lambda_nI_inverse(K, lambda_seq);
weights_seq = alpha_computation(KlambdaInv_seq, Vphi0);

end 

