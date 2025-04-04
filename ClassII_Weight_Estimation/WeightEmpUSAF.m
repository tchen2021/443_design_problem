function W_EmpennageUSAF = WeightEmpUSAF(Plane)
    %{
     Notes: 
     S is the area
     l is the distance from the wing quarter cord to the quarter
     cord of the horizontal/vertical stabilizer
     b is the tail span
     t_r is the max root thickness
    %}
    W_TO = Plane.W_dg;
    n_ult = Plane.Nz;
    
    S_h = Plane.S_h;
    l_h = Plane.l_h;
    b_h = Plane.b_h;
    t_r_h = Plane.t_r_h;
    S_v = Plane.S_v;
    l_v = Plane.l_v;
    b_v = Plane.b_v;
    t_r_v = Plane.t_r_v;

    W_h = 127 * (((W_TO * n_ult) / 10^5)^0.87 * (S_h / 100)^1.2 ...
          * 0.289 * (l_h / 10)^0.483 * (b_h / t_r_h)^(0.5))^0.458;

    W_v = 98.5 * (((W_TO * n_ult) / 10^5)^0.87 * (S_v / 100)^1.2 ...
          * 0.289 * (b_v / t_r_v)^(0.5))^0.483;

    W_EmpennageUSAF = W_h + W_v;
end