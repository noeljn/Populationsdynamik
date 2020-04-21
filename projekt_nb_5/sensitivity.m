function sensitivity(dataT2a, dataT2b)

    valuesV = [];
    func = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
    -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
    -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];

    %All values are changed to 120% the original value from func. One by one
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

    funcvec = {funca1,funca2,funca3,funcb1,funcb2,funcb3,funcc1,funcc2};    %Cell array with all functions
    
    %Calculating the last value for plants at t = 3 with Runge-Kutta 4
    %for every function where one coefficient from a1 - c2 is changed 
    %and storing the values for display together with the reference.
    reference = rk4(3,[dataT2a;2],dataT2b,func,2);
    for i = 1:8
        dataT3 = rk4(3,[dataT2a;2],dataT2b,funcvec{1,i},2);
        valuesV(end + 1) = dataT3(2,end);
    end

    disp("Running sensitivity of V");
    disp('Reference for sensitivity, no changes')
    disp(reference(2,end))
    disp("All values are changed to 120% the original value.");
    disp("Displayed one by one from a1 - c2")
    disp(valuesV');


end