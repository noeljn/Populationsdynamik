interpolT1()    %Finding T1

%Finding constant values for V and S and T2 after T1
fun = @(V,S)[15.*V+1.7*10^(-5).*V.^2-0.022.*V.*S;
             -1.9.*S.^(1.4)+0.088.*V.^(0.6).*S.^(0.8)];
jac = @(V,S)[15+3.4*10^(-5).*V-0.022.*S,-0.022.*V;
             0.0528.*S.^(0.8).*V.^(-0.4),-2.66.*S.^(0.4)+0.0704.*V.^(0.6).*S.^(-0.2)];
roots=newtonsys(10^(-7),[880000;700],fun,jac);
disp('V constant value:')
disp(roots(1))
disp('S constant value:')
disp(roots(2))