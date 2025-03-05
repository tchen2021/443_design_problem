
TUC = 6.7;
Wolv = 6.5;

[Tucano_e] = oswaldcalc(TUC)
[Wolver_e] = oswaldcalc(Wolv)


function [e] = oswaldcalc(plane)
e = 1.78*(1-(0.045*(plane^0.68)))-0.64;



end