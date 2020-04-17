function expansionMain3(dataT3)

%Golden ratio
r = (sqrt(5)-1)/2;
q = 1 - r;
%We can see on the graph from the plot with earlier that maximum is between
%0 and 0.1
a = 0;
b = 0.1;

x1 = a + q*(b-a);
data = harvesting(dataT3,x1,1);
index = find(data(1,:) == (7 + x1));
F1 = data(2,index - 1) - 100;

x2 = a + r*(b-a);
data = harvesting(dataT3,x2,1);
index = find(data(1,:) == (7 + x2));
F2 = data(2,index - 1) - 100;

count = 0;
while count < 100
    if F1 > F2
        b = x2;
        x2 = x1;
        F2 = F1;
        x1 = a + q*(b-a);
        data = harvesting(dataT3,x1,1);
        index = find(data(1,:) == (7 + x1));
        F1 = data(2,index - 1) - 100;
    else
        a = x1;
        x1 = x2;
        F1 = F2;
        x2 = a + r*(b-a);
        data = harvesting(dataT3,x2,1);
        index = find(data(1,:) == (7 + x2));
        F2 = data(2,index - 1) - 100;
    end
    count = count + 1;
end

fprintf("Maximum harvest is made when tau is %f\n",a);


%Calculates how much the harvest is year 7 depending on tau
collectedCrops = [];
for i = 0:0.01:0.99
    data = harvesting(dataT3, i, 1);    %With pestisides
    index = find(data(1,:) == (7 + i)); %Returns harvest for year 7
    harvest = data(2,index - 1) - 100;
    
    collectedCrops(1,end+1) = i;
    collectedCrops(2,end) = harvest;
    
    
    data = harvesting(dataT3, i, 0);    %Without pestisides
    index = find(data(1,:) == (7 + i));
    harvest = data(2,index - 1) - 100;
    collectedCrops(3,end) = harvest;
end

figure(5)
plot(collectedCrops(1,:),collectedCrops(2,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
plot(collectedCrops(1,:),collectedCrops(3,:),'color',[0.6666,0.2,0.3490]);
xlabel('Time')
ylabel('Harvest')
legend('With Pestaside','Without Pestaside')
title('Plot of harvest dependent on Tau year 7')


f1 = figure(10);
copyobj(get(figure(1),'children'),f1);
f2 = figure(11);
copyobj(get(figure(2),'children'),f2);
f3 = figure(12);
copyobj(get(figure(3),'children'),f3);
f4 = figure(13);
copyobj(get(figure(1),'children'),f4);
f5 = figure(14);
copyobj(get(figure(2),'children'),f5);
f6 = figure(15);
copyobj(get(figure(3),'children'),f6);

figurevec = [f1,f4,f2,f5,f3,f6];
tauvec = [0;0.2;0.5;0.8];
for n = 0:1
    for i = 1:4
        data = harvesting(dataT3,tauvec(i),n);
        figure(figurevec(1+n))
        plot(data(1,:),data(2,:));
        hold on
        figure(figurevec(3+n))
        plot(data(1,:),data(3,:));
        hold on
        figure(figurevec(5+n))
        plot(data(1,:),data(4,:));
        hold on
    end
    
    wop = "without pestisides";
    wp = "with pestisides";
    if n == 0
        title1 = strcat("Plants ",wop);
        title2 = strcat("Mice ",wop);
        title3 = strcat("Snakes ",wop);
    elseif n == 1
        title1 = strcat("Plants ",wp);
        title2 = strcat("Mice ",wp);
        title3 = strcat("Snakes ",wp);
    end
    figure(figurevec(1+n))
    xlabel('Time')
    ylabel('Population')
    title(title1)
    legend('Year 0-T1','Year T1-1.5','Year 1.5-3','Tau 0','Tau 0.2','Tau 0.5','Tau 0.8','Location','northeast')
    figure(figurevec(3+n))
    xlabel('Time')
    ylabel('Population')
    title(title2)
    legend('Year T1-1.5','Year 1.5-3','Tau 0','Tau 0.2','Tau 0.5','Tau 0.8','Location','northeast')
    figure(figurevec(5+n))
    xlabel('Time')
    ylabel('Population')
    title(title3)
    legend('Year 1.5-3','Tau 0','Tau 0.2','Tau 0.5','Tau 0.8','Location','northeast')
end

data = harvesting(dataT3,a,1);
figure(1)
plot(data(1,:),data(2,:),'color',[0.2980, 0, 0.6]);
hold on
figure(2)
plot(data(1,:),data(3,:),'color',[0.2980, 0, 0.6]);
hold on
figure(3)
plot(data(1,:),data(4,:),'color',[0.2980, 0, 0.6]);
hold on

end

function [data] = harvesting(dataT3,tau,pest)
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
    
    if tau ~= 0
        year = rk4(t + tau,start(2:4,end),start(1,end),func2,2); %Collecting data between t and tau
    else
        year = start;
    end
    
    %Harvest
    year(2,end) = 100;
    data = [data(:,1:end-1),year];
    start = [data(1:2,end);data(3,end);data(4,end)];
    
    
    %Pest
    t = t + 1;
    if pest == 1
        year = rk4(t,start(2:4,end),start(1,end),func2,2);
        data = [data(:,1:end-1),year];
        start = [data(1:2,end);data(3,end)*0.3;data(4,end)*0.8];
    elseif pest == 0
        year = rk4(t,start(2:4,end),start(1,end),func2,2);
        data = [data(:,1:end-1),year];
        start = [data(1:2,end);data(3,end);data(4,end)];
    end
    
end
end