function fvalue= spherical_cap_SD15_func(xyz)
    x = xyz(:, 1);
    y = xyz(:, 2);
    z = xyz(:, 3);
    center_points = SD(15);
    threshold = pi/32.0;
    num_centers = size(center_points,1);
    for i = 1:1:num_centers
         center_point = center_points(i,:);
        if i ==1
            fvalue = spherical_cap_single(x,y,z, center_point, threshold);
        else
            fvalue = fvalue + spherical_cap_single(x,y,z,center_point, threshold);
        end
    end
    fvalue = fvalue .* 1.0;
end

