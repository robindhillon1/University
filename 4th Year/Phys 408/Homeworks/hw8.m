lamb = 10.6E-4;

syms z w0
eqns = [(0.1699)^2 == w0^2*(1+(z/(pi*w0^2/lamb))^2), (0.3380)^2 == w0^2*(1+((z+10)/(pi*w0^2/lamb))^2)];
S = solve(eqns);

round(double(S.z), 10)

round(double(S.w0), 10)
