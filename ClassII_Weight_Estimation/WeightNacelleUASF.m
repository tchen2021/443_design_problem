function W_NacelleUSAF = WeightNacelleUASF(Plane)
    Plane.W_eng = W_eng; %Weight of engine group [lbs]
    Plane.Ne = Ne;       %Number of engines
    W_NacelleUSAF = 2.575*(W_eng^0.922)*Ne;
end

