function W_LandingGear = WeightLandingGear(plane)
    W_L = plane.W_dg*0.97;       %design landing weight in lbs see table 3.3 in part 1
    n_ult_1 = 5.7;   %ultimate load factor for landing
    l_s_m=plane.l_s_m;   %shock strut length of main gear super tucano
    %plane.l_s_n = l_s_n;   %shock strut length of nose gear

    % Compute W_g
    W_LandingGear = 0.054 * (l_s_m)^0.501 * (W_L * n_ult_1)^0.684;
end