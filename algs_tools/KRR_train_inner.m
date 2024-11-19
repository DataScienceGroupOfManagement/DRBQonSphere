function [alpha_seq, Ktr] =  KRR_train_inner(in_data, out_data, lambda_seq, KerPara)
% in_data - Input to the functio to be regressed.  N (points) x D (dimensional)
% out_data - Ouput of the function to be regressed. N x 1 (points)
% lambda_seq - Regularization Parameter Sequence

if nargin < 3
    lambda_seq = Lambda_q(1, 2, 10);
end
if nargin < 4
    KerPara.KernelType = 1;
end

if size(in_data,1) ~= size(out_data,1)
    fprintf('\nTotal number of points for function input and output are unequal');
    fprintf('\n Exitting program');
    return
else
    %% Compute K(x,x') on training set  
    Ktr = KernelComputation(in_data, in_data, KerPara);
    Ktr = (Ktr+Ktr')/2;
    %% Compute adaptive value of lambda
    KlambdaInv_seq = K_add_lambda_nI_inverse(Ktr, lambda_seq);
    alpha_seq = alpha_computation(KlambdaInv_seq, out_data);
end