function fvalue = small_func(xyz,left_threshold, right_thrshold)
x = xyz(:, 1);
y = xyz(:, 2);
z = xyz(:, 3);

theta = acos(z);
phi = atan(y./x);

temp_index = theta >=left_threshold  & theta <= right_thrshold ;
% temp_index
fvalue = temp_index .* 1.0;

end
