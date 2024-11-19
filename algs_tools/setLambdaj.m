function V = setLambdaj(D, c0)

N = size(D, 1);
idx_rand = randperm(N);
idx_tmp = idx_rand(1);
V = zeros(1, N);
V(1) = idx_tmp;
Vl = find(D(idx_tmp,:) > c0);

count = 1;
while ~isempty(Vl)
    count = count + 1;
    D_tmp = D(Vl, Vl);
    N_tmp = size(D_tmp, 1);
    idx_rand = randperm(N_tmp);
    idx_tmp = idx_rand(1);
    V_tmp = Vl(idx_tmp);
    V(count) = V_tmp;

    Vl_tmp = D_tmp(idx_tmp,:) > c0;
    Vl = Vl(Vl_tmp);   
end

V(V==0) = [];