clear
format long

%Finding T1 and plotting the curve until T1
vMax = 15/(1.7*10^(-5));        %Finding the constant value
fun = @(v) 15*v - 1.7*10^(-5)*v^2;
v95 = 0.95*vMax;
y = rk4(vMax*0.99,100,0,fun,1); %Retriving plant population with function 
                                %for Runge-Kutta 4
[dataT1,T1] = interpolT1(y,v95);
disp('T1 is:')
disp(T1)
figure(1)
plot(dataT1(1,:),dataT1(2,:),'-b')
hold on

%Finding constant values for V and S
func1 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)];
jac1 = @(u)[15-3.4*10^(-5).*u(1)-0.022.*u(2),-0.022.*u(1);
             0.0528.*u(2).^(0.8).*u(1).^(-0.4),-2.66.*u(2).^(0.4)+0.0704.*u(1).^(0.6).*u(2).^(-0.2)];
constants=newtonsys([100000;700],func1,jac1);

%Collecting data with Runge-Kutta 4 and 
%Plotting the population growth T1-T2
dataT2 = rk4(1.5,[v95;2],T1,func1,2);

figure(1)
plot(dataT2(1,:),dataT2(2,:),'color',[0.6666,0.2,0.3490])
hold on
figure(2)
plot(dataT2(1,:),dataT2(3,:),'color',[0.6666,0.2,0.3490])
hold on


%Calculating the difference in permille 
diffV = abs(dataT2(2,end) - constants(1))/constants(1)*1000;
diffS = abs(dataT2(3,end) - constants(2))/constants(2)*1000;
Tab1 = table(constants(1),constants(2),diffV,diffS);
Tab1.Properties.VariableNames={'V_constant','S_constant','Diff_V_permille','Diff_S_permille'};
disp(Tab1)


%Finding new constant values for V, S and R
func2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
jac2 = @(u)[15-3.4*10^(-5).*u(1)-0.022.*u(2),-0.022.*u(1),0;
         0.0528.*u(2).^(0.8).*u(1).^(-0.4),-2.66.*u(2).^(0.4)+0.0704.*u(1).^(0.6).*u(2).^(-0.2)-1.4*u(3),-1.4.*u(2);
         0,0.028.*sqrt(u(3)),-1.8 + 0.028.*u(2).*0.5.*1./sqrt(u(3))];

constants = newtonsys([100000;600;20],func2,jac2);

%Collecting data between T2-T3 with Runge-Kutta 4
dataT3 = rk4(3,[dataT2(2:3,end);2],dataT2(1,end),func2,2);

%Plotting the population growth T2-T3

figure(1)
plot(dataT3(1,:),dataT3(2,:),'color',[0,0.8784,0.5294])
hold on
figure(2)
plot(dataT3(1,:),dataT3(3,:),'color',[0,0.8784,0.5294])
hold on
figure(3)
plot(dataT3(1,:),dataT3(4,:),'color',[0,0.8784,0.5294])
hold on


%Tab2 with calculated constant when T=inf
Tab2 = table(constants(1),constants(2),constants(3));
Tab2.Properties.VariableNames={'V_constant','S_constant','R_constant'};
disp(Tab2)
%Tab3 with calculated data Runge-Kutta 4
Tab3 = table(dataT3(2,end),dataT3(3,end),dataT3(4,end));
Tab3.Properties.VariableNames={'V_year_3','S_year_3','R_year_3'};
disp(Tab3);

%Calculating the effects of spraying and plotting the results
start = [dataT3(1:2,end);dataT3(3,end)*0.3;dataT3(4,end)*0.8];
data = start;
for i=1:5
    year = rk4(3+i,start(2:4,end),start(1,end),func2,2);
    data = [data(:,1:end-1),year];
    start = [data(1:2,end);data(3,end)*0.3;data(4,end)*0.8];
end

figure(1)
plot(data(1,:),data(2,:),'color',[1, 0.5490, 0.5176]);
hold on
figure(2)
plot(data(1,:),data(3,:),'color',[1, 0.5490, 0.5176]);
hold on
figure(3)
plot(data(1,:),data(4,:),'color',[1, 0.5490, 0.5176]);
hold on

%Setting labels on the figures

figure(1)
xlabel('Time')
ylabel('Population')
legend('0-T1','T1-T2','T2-T3','Spraying','Location','northeast')
title('Plants')

figure(2)
xlabel('Time')
ylabel('Population')
legend('T1-T2','T2-T3','Spraying','Location','northeast')
title('Mice')

figure(3)
xlabel('Time')
ylabel('Population')
legend('T2-T3','Spraying','Location','southeast')
title('Snakes')

%Sensitivity
%sensitivity(dataT2(2:3,end),dataT2(1,end));




%{
%Expansion
expansionMain(dataT3);
%reliabilityExpansion
reliabilityExpansion(dataT3);
%}







