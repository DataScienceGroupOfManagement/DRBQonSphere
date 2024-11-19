function fvalue= sim_func_circles(xyz)
    interval = pi/16;
    fvalue = small_func(xyz, 0, 0 + interval)  + small_func(xyz, pi/8, pi/8+ interval) + small_func(xyz, pi/4, pi/4+ interval);
end

