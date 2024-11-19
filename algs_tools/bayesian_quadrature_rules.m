function [integrate_value_seq ] = bayesian_quadrature_rules(out_data, weights_seq)
    % out_data - Ouput of the function to be integrated. N x 1 (points)
    % weights_seq - bayesian quadrature weights. N x length(lambda_seq) 
    integrate_value_seq = out_data' * weights_seq; 
end 

