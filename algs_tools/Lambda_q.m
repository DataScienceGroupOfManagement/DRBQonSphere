function lambda_q_seq = Lambda_q(kappa2, q, m)
q_seq = q.^(1:m);
lambda_q_seq = kappa2 ./ q_seq;