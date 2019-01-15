%x_min should be higher than 0 
%x_max lower than 110
%ax= 0.5 for x = 1-109

%if UseAccumulatedEndGroup =1 then for x_max a accumulated group x_max+ was
%used in the data Matrix M=m(x,t)
%for q(x,t) there was no adjustment
%which means last row of M is not comparable with q in this case
if UseAccumulatedEndGroup==1
     disp('Beware, the accumulated end age used is not comparable with Cairns')
end

%Transfrom q_C estimates by Cairns fit to estimates for death rates M

%from HMD: qx=mx/(1+(1-ax)*mx), with ax=0.5 for x=1-109, WRONG said STADJE
M_HMD= q_C./(1-(1-0.5)*q_C);
%H= abs(M-q./(1-(1-0.5)*q)); %difference Matrix only from Transforming HMD data
%sum(sum(H)) %total absolut difference, should be 0?

%under piecewise constant assumption (2.1): qx=1-exp(-mx),mx=-log(1-qx),
%TRUE said stadtje
M_Cairns=-log(1-q_C);
%H_21=abs(-log(1-q)-M); %difference Matrix only from Transforming HMD data
%sum(sum(H_21))%total absolut difference 

figure('name','m(x,t) Fits of Lee-Carter Modell and Cairns Modell');
ax1 = subplot(2,1,1);
%plots development over time for a specific age
age=100; %beetwen x_min_c and x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);
time_C=linspace(t_min_C,t_max_C,t_C);
ages_C=linspace(x_min_C,x_max_C,x_C);
time_full=linspace(t_start,t_end,t_full);
ages_full=linspace(0,x_full-1,x_full);

%x_min in zeile 1
%x_min +1 in zeile 2
%x_min +y=age in zeile y+1, y=age-x_min
y=age-x_min+1;
y_C=age-x_min_C+1;

tempM=M_full(age+1,:);
index42=tempM~=42;

plot(time_full(index42),tempM(index42),'g')
hold on
plot(time,M_LC(y,:),'r')
hold on
plot(time,M_SVD(y,:),'k')
hold on
plot(time,M_ML(y,:),'b')
hold on
plot(time,M_wLC(y,:),'c')
hold on
plot(time,M_wLC2(y,:),'m')
hold on
plot(time_C,M_Cairns(y_C,:),'y')

h=num2str(age);
title(ax1,['m(x,t) over time for age ' h])
legend('historical','LC','SVD','ML','wLC','wLC2','Cairns')

ax2 = subplot(2,1,2);

%plots development over age for a specific year

%for this part time range for Cairns and LC is assumed to be same
%z.B. 1956-2015

year=t; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
%view_start=1; %1 for start age x_min %%%%%%%%%%%%%%%%%%%
%view_end=x; %x for end age x-1, i.e. 100+ for x=101 %%%%%%%%%%%
%view=linspace(view_start,view_end,view_end-view_start+1);

%lnM_full=log(M_full); contains 0 and -1 
lnM_Cairns=log(M_Cairns);

templnM=lnM(:,year);
indexln42=templnM~=log(42);

plot(ages(indexln42),templnM(indexln42),'g')
hold on
plot(ages,lnM_LC(:,year),'r')
hold on
plot(ages,lnM_SVD(:,year),'k')
hold on
plot(ages,lnM_ML(:,year),'b')
hold on
plot(ages,lnM_wLC(:,year),'c')
hold on
plot(ages,lnM_wLC2(:,year),'m')
hold on
plot(ages_C,lnM_Cairns(:,year),'y')
h=num2str(year+t_min-1);
title(ax2,['ln[m(x,t)] over ages in year ' h])

%clearvars h h1 h2 n age ages ax1 ax2 i j time time_C view view_end view_start
%clearvars year y y_C time_full ages_C ages_full tempM index42 templnM
%clearvars indexln42