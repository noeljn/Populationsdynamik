clear
format long

%Finding T1
T1 = interpolT1();

    
%Finding constant values for V and S
fun = @(V,S)[15.*V-1.7*10^(-5).*V.^2-0.022.*V.*S;
             -1.9.*S.^(1.4)+0.088.*V.^(0.6).*S.^(0.8)];
jac = @(V,S)[15-3.4*10^(-5).*V-0.022.*S,-0.022.*V;
             0.0528.*S.^(0.8).*V.^(-0.4),-2.66.*S.^(0.4)+0.0704.*V.^(0.6).*S.^(-0.2)];
constants=newtonsys(10^(-6),[100000;700],fun,jac);


T2 = punkt4();
diffV = (T2(1) - constants(1))/constants(1)*1000;
diffS = (T2(2) - constants(2))/constants(2)*1000;
T = table(T1,constants(1),constants(2),diffV,diffS);
T.Properties.VariableNames={'T1','V_constant','S_constant','Diff_V‰','Diff_S‰'};


disp(T);
