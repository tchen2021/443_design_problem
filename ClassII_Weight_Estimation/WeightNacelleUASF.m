function W_NacelleUSAF = WeightNacelleUASF(Plane)
    W_eng = Plane.W_eng; %Weight of engine group [lbs]
    Ne = 1;       %Number of engines
    W_NacelleUSAF = 2.575*(W_eng^0.922)*Ne;
end

