function D = geodesic_distance_compute(X, Y)
% https://math.stackexchange.com/questions/1304169/distance-between-two-points-on-a-sphere
Z = X * Y';
D = acos(Z);
D = real(D);
D = D-diag(diag(D));
