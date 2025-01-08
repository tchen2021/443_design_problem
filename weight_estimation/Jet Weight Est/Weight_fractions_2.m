function [W_TO, W_E, W_F,MFF] = Weight_fractions_2(c_j, W, W_pl, Range, E, V, L_D, EWF, RFF, initial_guess,p)
    %future adjustment: have function take in parameters as vectors, 1
    %value for cruise and 1 for loiter
    %assumes same L_D during loiter and cruise

    %take velocity as array, first number is cruise, second is loiter,
    %add different conditions
    %Range(1) = cruise, Range(2) = attac
    if p == 1
        W(5) = (exp(Range(1)/( (V(1)/c_j)*(L_D(1)) )))^-1;
        if E == 0
            W(6) = 1;
        else
            W(6) =  (1 / ( exp(( (E+0.5)  * c_j )/L_D(2) ) ) );
        end
        W(7) = (exp(Range(1)/( (V(1)/c_j)*(L_D(1)) )))^-1;
        MFF = prod(W,"all");
        
        W_TO0 = initial_guess(1);
        W_TO = W_TO0;
        W_E_max = EWF*W_TO;
     
        fun  = @(W_TO) (W_TO - W_TO*(1-MFF)- W_pl) - W_E_max;
        W_TO = fzero(fun, W_TO0);
        W_F = (1-MFF)*W_TO;
        W_E = EWF * W_TO;
        
    elseif p == 2
        W(5) = (exp(Range(1)/( (V(2)/c_j)*(L_D(1)) )))^-1;
        if E == 0
            W(6) = 1;
        else
            W(6) =  (1 / ( exp(( (E+0.5)  * c_j )/L_D(2) ) ) );
        end

        W(8) = (exp(Range(2)/( (V(2)/c_j)*(L_D(1)) )))^-1;
        W(10) = (exp(Range(1)/( (V(1)/c_j)*(L_D(1)) )))^-1;
             
        MFF = prod(W,"all");
        
        W_TO0 = initial_guess(1);
        W_TO = W_TO0;
        W_E_max = EWF*W_TO;
      
        fun  = @(W_TO) (W_TO - W_TO*(1-MFF)- W_pl) - W_E_max;
        W_TO = fzero(fun, W_TO0);
        W_F = (1-MFF)*W_TO;
        W_E = EWF * W_TO;
    end   
end
