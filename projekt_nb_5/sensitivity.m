function sensitivity(dataT2a, dataT2b)

valuesV = [];
  %All values are changed to 120% the original value. One by one
         funca1 = @(u)[18.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funca2 = @(u)[15.*u(1)-2.04*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funca3 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.0264.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funcb1 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -2.28.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funcb2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.1056.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funcb3 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.68.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funcc1 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -2.16.*u(3) + 0.028.*u(2).*sqrt(u(3))];
         funcc2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.0336.*u(2).*sqrt(u(3))];

dataT3 = rk4(3,[dataT2a;2],dataT2b,funca1,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funca2,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funca3,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funcb1,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funcb2,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funcb3,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funcc1,2);
valuesV(end + 1) = dataT3(2,end);

dataT3 = rk4(3,[dataT2a;2],dataT2b,funcc2,2);
valuesV(end + 1) = dataT3(2,end);
disp("Sensitivity of V");
disp("All values are changed to 120% the original value. One by one form a1 - c2");
disp(valuesV');


end