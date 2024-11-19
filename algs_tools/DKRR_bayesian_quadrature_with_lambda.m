function [quadrature_value, mean_time] = DKRR_bayesian_quadrature_with_lambda(train_x,train_y, KerPara, lambda_value, m,  phi0_value)
    
    Machine = cell(1,m);
    lambda_optimal = lambda_value;
    
    Ntotal = size(train_x,1);
    % Data division
    Nvec = zeros(m,1);
    n = floor(Ntotal/m);
    for k = 1:m-1
        Nvec(k) = n;
    end
    Nvec(m) = Ntotal-n*(m-1);
    indjCell =  cell(1,m);
    count = 0;
    for k = 1:m
        indjCell{k} = count+1:count+Nvec(k);
        count = count+Nvec(k);
    end
    data_idx_shuffle = randperm(Ntotal); 
    for k = 1:m
        Machine{k}.IndData = data_idx_shuffle(indjCell{k});
        Machine{k}.train_x = train_x(Machine{k}.IndData,:);
        Machine{k}.train_y = train_y(Machine{k}.IndData,:);
        Machine{k}.num_samples = size(Machine{k}.train_x,1);
        Machine{k}.Vphi0 = ones(Machine{k}.num_samples,1)*phi0_value;
    end
    
    % local processing
    t_local_time = zeros(m,1);
    for k = 1:m
        tic;
        Machine{k}.baye_weights =  KRR_train_inner(Machine{k}.train_x, Machine{k}.Vphi0, lambda_optimal, KerPara);
        t_tmp = toc;
        t_local_time(k) = t_tmp;
    end
    mean_time = mean(t_local_time);
    
    % synthesization
    quadrature_value = 0;
    for k = 1:m
        quadrature_value = quadrature_value + (Machine{k}.num_samples / Ntotal) * bayesian_quadrature_rules(Machine{k}.train_y, Machine{k}.baye_weights);
    end
end