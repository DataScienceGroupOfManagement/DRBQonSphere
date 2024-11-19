function Vs = setsLambda(D, c0)

V_left = 1:size(D,1);
D_tmp = D;
Vs = cell(1,1);

count = 0;
while ~isempty(V_left)
    count = count + 1;
    V_tmp = setLambdaj(D_tmp, c0);
    V = V_left(V_tmp);
    Vs{count} = V;

    V_left = setdiff(V_left, V);
    D_tmp = D(V_left, V_left);
end

