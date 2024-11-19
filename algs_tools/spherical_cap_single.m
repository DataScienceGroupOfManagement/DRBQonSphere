function  fvalue= spherical_cap_single(x, y, z, center_point, threshold)
    s1s2 = sqrt((x-center_point(1)).^2 +  (y-center_point(2)).^2 + (z-center_point(3)).^2);
    r = acos(1-0.5.*s1s2.^2);
    temp_index = r<=threshold;
    fvalue = temp_index .* 1.0;
end