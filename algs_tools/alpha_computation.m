function alpha_seq = alpha_computation(KlambdaInv_seq, y)
lambda_num = size(KlambdaInv_seq, 3);
N = length(y);
alpha_seq = zeros(N, lambda_num);
for k = 1:lambda_num
    alpha_seq(:, k) = KlambdaInv_seq(:,:,k)*y;
end