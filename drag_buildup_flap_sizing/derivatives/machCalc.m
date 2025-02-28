function M = machCalc(V,h)

V = V* 1.68781;     %convert from knots to ft/s
[~, ~, ~, a] = atmosphere(h);   %obtain speed of sound
M = V/a;

end