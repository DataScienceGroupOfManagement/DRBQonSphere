function [Nvec, Vs_re] = setsReLambda(D, Vs, c0)
N = size(D, 1);
N_ave = floor(N/length(Vs));

Nvec_bef = zeros(1, length(Vs));
for i = 1:length(Vs)
    Nvec_bef(i) = length(Vs{i});
end

idx_need_add = find(Nvec_bef < N_ave);
idx_need_minus = find(Nvec_bef > N_ave);

V_pre = cell(1, length(idx_need_minus));
V_adj = cell(1, length(idx_need_minus));
for i = 1:length(idx_need_minus)
    V_tmp = Vs{idx_need_minus(i)};
    N_tmp = Nvec_bef(idx_need_minus(i));
    
    idx_rand = randperm(N_tmp);
    V_pre{i} = V_tmp(idx_rand(1:N_ave));
    V_adj{i} = setdiff(V_tmp, V_pre{i});
end

Vs_re = Vs;
for i = length(idx_need_add):-1:1
    idx_set_need_add = idx_need_add(i);
    idx_tmp = Vs{idx_set_need_add};
    
    for j = 1:length(idx_need_minus)
%         idx_adj_tmp = V_adj{idx_need_minus(j)};
        idx_adj_tmp = V_adj{j};
        if ~isempty(idx_adj_tmp)
            D_tmp = D(idx_tmp, idx_adj_tmp);
            idx_add_tmp1 = sum(D_tmp>c0, 1) == length(idx_tmp);
            idx_add_tmp = idx_adj_tmp(idx_add_tmp1);
            N_need_add = N_ave - length(idx_tmp); 
            if length(idx_add_tmp) >= N_need_add
                idx_tmp = [idx_tmp idx_add_tmp(1:N_need_add)];
                V_adj{j} = setdiff(idx_adj_tmp, idx_add_tmp(1:N_need_add));
                Vs_re{idx_set_need_add} = idx_tmp;
                break;
            else
                idx_tmp = [idx_tmp idx_add_tmp];
                V_adj{j} = setdiff(idx_adj_tmp, idx_add_tmp);
            end
        end
    end
    Vs_re{idx_set_need_add} = idx_tmp;
end

for i = 1:length(idx_need_minus)
    Vs_re{idx_need_minus(i)} = [V_pre{i} V_adj{i}];
end

Nvec = zeros(1, length(Vs_re));
for i = 1:length(Vs_re)
    Nvec(i) = length(Vs_re{i});
end



