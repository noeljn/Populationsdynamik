function expansionMain(dataT3)

%Golden ratio
r = (sqrt(5)-1)/2;
q = 1 - r;
%We can se on the graph from the plot with earlier that maximum is between
%0 and 0.1
a = 0;
b = 0.1;

x1 = a + q*(b-a);
data = expansion(dataT3,x1,1);
index = find(data(1,:) == (7 + x1));
F1 = data(2,index - 1) - 100;

x2 = a + r*(b-a);
data = expansion(dataT3,x2,1);
index = find(data(1,:) == (7 + x2));
F2 = data(2,index - 1) - 100;

count = 0;
while count < 100
    if F1 > F2
        b = x2;
        x2 = x1;
        F2 = F1;
        x1 = a + q*(b-a);
        data = expansion(dataT3,x1,1);
        index = find(data(1,:) == (7 + x1));
        F1 = data(2,index - 1) - 100;
    else
        a = x1;
        x1 = x2;
        F1 = F2;
        x2 = a + r*(b-a);
        data = expansion(dataT3,x2,1);
        index = find(data(1,:) == (7 + x2));
        F2 = data(2,index - 1) - 100;
    end
    count = count + 1;
end
Tab = table(a);
Tab.Properties.VariableNames={'Max_Harvest(tau)'};
disp(Tab)


%{
%Plot V S T for t ∈ [0, 8] with and without pestaside
for i = 0:1
    %tau = 0
    data = expansion(dataT3, 0, i);
    figure(1)
    plot(data(1,:),data(2,:),'color',[0.9568, 0.5019, 0.1254]);
    hold on
    figure(2)
    plot(data(1,:),data(3,:),'color',[0.9568, 0.5019, 0.1254]);
    hold on
    figure(3)
    plot(data(1,:),data(4,:),'color',[0.9568, 0.5019, 0.1254]);
    hold on

    %tau = 0.5
    data = expansion(dataT3, 0.5, i);
    figure(1)
    plot(data(1,:),data(2,:),'color',[0.9568, 0.5019, 0.1254]);
    hold on
    figure(2)
    plot(data(1,:),data(3,:),'color',[0.9568, 0.5019, 0.1254]);
    hold on
    figure(3)
    plot(data(1,:),data(4,:),'color',[0.9568, 0.5019, 0.1254]);
    hold on
end
%}
%Collected Crops over time
collectedCrops = [];
for i = 0.01:0.01:0.99
    data = expansion(dataT3, i, 1);
    %Returns harvest for year 7
    index = find(data(1,:) == (7 + i));
    harvest = data(2,index - 1) - 100;
    
    collectedCrops(2,end + 1) = harvest;
    collectedCrops(1,end) = i;
    
    data = expansion(dataT3, i, 0);
    index = find(data(1,:) == (7 + i));
    harvest = data(2,index - 1) - 100;
    collectedCrops(3,end) = harvest;
end

figure(4)
plot(collectedCrops(1,:),collectedCrops(2,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
plot(collectedCrops(1,:),collectedCrops(3,:),'color',[0.6666,0.2,0.3490]);

figure(4)
xlabel('Time')
ylabel('Harvest')
legend('With Pestaside','Without Pestaside')
title('Plot of harvest dependent on Tau')

end


function [data] = expansion(dataT3,tau,pest)
func2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];

%Collecting data between T3-T8 with Runge-Kutta 4 
t = 3;
if pest == 1
    start = [dataT3(1:2,end);dataT3(3,end)*0.3;dataT3(4,end)*0.8];
else
    start = [dataT3(1:2,end);dataT3(3,end);dataT3(4,end)];
end
data = start;

while t < 8
    
    %Collecting data between t and tau
    if tau ~= 0
        year = rk4(t + tau,start(2:4,end),start(1,end),func2,2);
        %Harvest
        year(2,end) = 100;
        data = [data(:,1:end-1),year];
        start = [data(1:2,end);data(3,end);data(4,end)];
    end
    
    %Pest
    t = t + 1;
    if pest == 1
        year = rk4(t,start(2:4,end),start(1,end),func2,2);
        data = [data(:,1:end-1),year];
        start = [data(1:2,end);data(3,end)*0.3;data(4,end)*0.8];
    end
    
end


%{
figure(1)
plot(data(1,:),data(2,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
figure(2)
plot(data(1,:),data(3,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
figure(3)
plot(data(1,:),data(4,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
%}


end