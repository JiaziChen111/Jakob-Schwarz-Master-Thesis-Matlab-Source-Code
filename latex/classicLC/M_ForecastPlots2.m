clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=110; %new maximum age group 100+, set to prevent 0 and NaN in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for fit
t_max=1996; %last year of data included for fit
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %do not replace 0 and -1 with 42 in M, instead here we set x_max=104

%load data with above parameters
BuildCentralDeathMatrix; 

W=zeros(x,t);
for i=0:x_max
    for j=1:t
    W(i+1,j)=t_min+j-1+x_max-i;
    end
end
 B=lnM==log(42);
    W(B)=0;


[a_W,b_W,k_W,steps_W,fval_W]=Fit_wLC(lnM,W);
[a_W2,b_W2,k_W2]=LeeCarterFit(D,E,a_W,b_W,k_W);

[a_wLC,b_wLC,k_wLC,steps_wLC,fval_wLC]=Fit_wLC(lnM);
[a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);

[a_D,b_D,k_D,steps_D,fval_D]=Fit_wLC(lnM,D);
[a_D2,b_D2,k_D2]=LeeCarterFit(D,E,a_D,b_D,k_D);

[a_ML,b_ML,k_ML,fval_ML]=PoissonFit(D, E);
[a_ML2,b_ML2,k_ML2]=LeeCarterFit(D,E,a_ML,b_ML,k_ML);
        
M_wLC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC(i,j)=exp(a_wLC(i)+b_wLC(i)*k_wLC(j));
        end
    end
M_wLC2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC2(i,j)=exp(a_wLC2(i)+b_wLC2(i)*k_wLC2(j));
        end
    end    
M_D=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D(i,j)=exp(a_D(i)+b_D(i)*k_D(j));
        end
    end    
    M_D2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D2(i,j)=exp(a_D2(i)+b_D2(i)*k_D2(j));
        end
    end    
  M_W=zeros(x,t);
    for i=1:x
        for j=1:t
            M_W(i,j)=exp(a_W(i)+b_W(i)*k_W(j));
        end
    end
M_W2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_W2(i,j)=exp(a_W2(i)+b_W2(i)*k_W2(j));
        end
    end
M_ML=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML(i,j)=exp(a_ML(i)+b_ML(i)*k_ML(j));
        end
    end
M_ML2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML2(i,j)=exp(a_ML2(i)+b_ML2(i)*k_ML2(j));
        end
    end
    Arima010=arima(0,1,0);
    Arima011=arima(0,1,1);
    Arima110=arima(1,1,0);
    Arima111=arima(1,1,1);
    Arima012=arima(0,1,2);
    Arima112=arima(1,1,2);
    Arima210=arima(2,1,0);
    Arima211=arima(2,1,1);
    Arima212=arima(2,1,2);
   
     t_future=26;
     n_paths=3;
    
 [Est_wLC,~,res_wLC,Mean_wLC,UB_wLC,LB_wLC, Ysim_wLC ,simU_wLC, simL_wLC, YsimParUnc_wLC, U_ParUnc_wLC, L_ParUnc_wLC] = ArimaForecast(t_future,n_paths,Arima011,k_wLC);  
 [Est_wLC2,~,res_wLC2,Mean_wLC2,UB_wLC2,LB_wLC2, Ysim_wLC2 ,simU_wLC2, simL_wLC2, YsimParUnc_wLC2, U_ParUnc_wLC2, L_ParUnc_wLC2] = ArimaForecast(t_future,n_paths,Arima011,k_wLC2);  
 [Est_W,EstCov_W,res_W,Mean_W,UB_W,LB_W, Ysim_W ,simU_W, simL_W, YsimParUnc_W, U_ParUnc_W, L_ParUnc_W] = ArimaForecast(t_future,n_paths,Arima010,k_W);  
 [Est_W2,~,res_W2,Mean_W2,UB_W2,LB_W2, Ysim_W2 ,simU_W2, simL_W2, YsimParUnc_W2, U_ParUnc_W2, L_ParUnc_W2] = ArimaForecast(t_future,n_paths,Arima011,k_W2);  
 [Est_ML,~,res_ML,Mean_ML,UB_ML,LB_ML, Ysim_ML ,simU_ML, simL_ML, YsimParUnc_ML, U_ParUnc_ML, L_ParUnc_ML] = ArimaForecast(t_future,n_paths,Arima011,k_ML);  
 [Est_ML2,~,res_ML2,Mean_ML2,UB_ML2,LB_ML2, Ysim_ML2 ,simU_ML2, simL_ML2, YsimParUnc_ML2, U_ParUnc_ML2, L_ParUnc_ML2] = ArimaForecast(t_future,n_paths,Arima011,k_ML2);  
 [Est_D,~,res_D,Mean_D,UB_D,LB_D, Ysim_D ,simU_D, simL_D, YsimParUnc_D, U_ParUnc_D, L_ParUnc_D] = ArimaForecast(t_future,n_paths,Arima011,k_D);  
 [Est_D2,~,res_D2,Mean_D2,UB_D2,LB_D2, Ysim_D2 ,simU_D2, simL_D2, YsimParUnc_D2, U_ParUnc_D2, L_ParUnc_D2] = ArimaForecast(t_future,n_paths,Arima011,k_D2);  
 
M_wLC_Forecast=zeros(x,t_future);
M_wLC_Forecast_UB=zeros(x,t_future);
M_wLC_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_wLC_Forecast(i,j)=exp(a_wLC(i)+b_wLC(i)*Mean_wLC(j));
             M_wLC_Forecast_UB(i,j)=exp(a_wLC(i)+b_wLC(i)*UB_wLC(j));
              M_wLC_Forecast_LB(i,j)=exp(a_wLC(i)+b_wLC(i)*LB_wLC(j));
        end
    end

M_ML_Forecast=zeros(x,t_future);
M_ML_Forecast_UB=zeros(x,t_future);
M_ML_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_ML_Forecast(i,j)=exp(a_ML(i)+b_ML(i)*Mean_ML(j));
             M_ML_Forecast_UB(i,j)=exp(a_ML(i)+b_ML(i)*UB_ML(j));
              M_ML_Forecast_LB(i,j)=exp(a_ML(i)+b_ML(i)*LB_ML(j));
        end
    end
   
M_ML2_Forecast=zeros(x,t_future);
M_ML2_Forecast_UB=zeros(x,t_future);
M_ML2_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_ML2_Forecast(i,j)=exp(a_ML2(i)+b_ML2(i)*Mean_ML2(j));
             M_ML2_Forecast_UB(i,j)=exp(a_ML2(i)+b_ML2(i)*UB_ML2(j));
              M_ML2_Forecast_LB(i,j)=exp(a_ML2(i)+b_ML2(i)*LB_ML2(j));
        end
    end
    
M_D_Forecast=zeros(x,t_future);
M_D_Forecast_UB=zeros(x,t_future);
M_D_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_D_Forecast(i,j)=exp(a_D(i)+b_D(i)*Mean_D(j));
             M_D_Forecast_UB(i,j)=exp(a_D(i)+b_D(i)*UB_D(j));
              M_D_Forecast_LB(i,j)=exp(a_D(i)+b_D(i)*LB_D(j));
        end
    end
 M_D2_Forecast=zeros(x,t_future);
M_D2_Forecast_UB=zeros(x,t_future);
M_D2_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_D2_Forecast(i,j)=exp(a_D2(i)+b_D2(i)*Mean_D2(j));
             M_D2_Forecast_UB(i,j)=exp(a_D2(i)+b_D2(i)*UB_D2(j));
              M_D2_Forecast_LB(i,j)=exp(a_D2(i)+b_D2(i)*LB_D2(j));
        end
    end
   M_W2_Forecast=zeros(x,t_future);
M_W2_Forecast_UB=zeros(x,t_future);
M_W2_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_W2_Forecast(i,j)=exp(a_W2(i)+b_W2(i)*Mean_W2(j));
             M_W2_Forecast_UB(i,j)=exp(a_W2(i)+b_W2(i)*UB_W2(j));
              M_W2_Forecast_LB(i,j)=exp(a_W2(i)+b_W2(i)*LB_W2(j));
        end
    end
     M_W_Forecast=zeros(x,t_future);
M_W_Forecast_UB=zeros(x,t_future);
M_W_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_W_Forecast(i,j)=exp(a_W(i)+b_W(i)*Mean_W(j));
             M_W_Forecast_UB(i,j)=exp(a_W(i)+b_W(i)*UB_W(j));
              M_W_Forecast_LB(i,j)=exp(a_W(i)+b_W(i)*LB_W(j));
        end
    end
      
M_wLC2_Forecast=zeros(x,t_future);
M_wLC2_Forecast_UB=zeros(x,t_future);
M_wLC2_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_wLC2_Forecast(i,j)=exp(a_wLC2(i)+b_wLC2(i)*Mean_wLC2(j));
             M_wLC2_Forecast_UB(i,j)=exp(a_wLC2(i)+b_wLC2(i)*UB_wLC2(j));
              M_wLC2_Forecast_LB(i,j)=exp(a_wLC2(i)+b_wLC2(i)*LB_wLC2(j));
        end
    end

lnM_W_Forecast=log(M_W_Forecast);
lnM_W_Forecast_UB=log(M_W_Forecast);
lnM_W_Forecast_LB=log(M_W_Forecast);

lnM_W2_Forecast=log(M_W2_Forecast);
lnM_W2_Forecast_UB=log(M_W2_Forecast);
lnM_W2_Forecast_LB=log(M_W2_Forecast);

lnM_ML_Forecast=log(M_ML_Forecast);
lnM_ML_Forecast_UB=log(M_ML_Forecast_UB);
lnM_ML_Forecast_LB=log(M_ML_Forecast_LB);

lnM_ML2_Forecast=log(M_ML2_Forecast);
lnM_ML2_Forecast_UB=log(M_ML2_Forecast_UB);
lnM_ML2_Forecast_LB=log(M_ML2_Forecast_LB);

lnM_wLC_Forecast=log(M_wLC_Forecast);
lnM_wLC_Forecast_UB=log(M_wLC_Forecast_UB);
lnM_wLC_Forecast_LB=log(M_wLC_Forecast_LB);

lnM_wLC2_Forecast=log(M_wLC2_Forecast);
lnM_wLC2_Forecast_UB=log(M_wLC2_Forecast_UB);
lnM_wLC2_Forecast_LB=log(M_wLC2_Forecast_LB);

lnM_D_Forecast=log(M_D_Forecast);
lnM_D_Forecast_UB=log(M_D_Forecast_UB);
lnM_D_Forecast_LB=log(M_D_Forecast_LB);

lnM_D2_Forecast=log(M_D2_Forecast);
lnM_D2_Forecast_UB=log(M_D2_Forecast_UB);
lnM_D2_Forecast_LB=log(M_D2_Forecast_LB);

lnM_full=log(M_full);    
    
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);   
time_future=linspace(t_max+1,t_max+t_future,t_future);

forecast1=figure('name','m(x,t) Forecast of Lee-Carter Modell');
ax1 = subplot(4,1,1);
%plots development over time for a specific age
age=0; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose

tempM=M_full(age+1,:);
index42=tempM~=42;
time_full=linspace(t_start,t_end,t_full);

plot(time_full(index42),tempM(index42),'g')
%hold on
%plot(time,M_wLC(age+1,:),'r')
%hold on
%plot(time,M_D(age+1,:),'k')
%hold on
%plot(time,M_ML(age+1,:),'b')
%hold on
%plot(time,M_W(age+1,:),'c')
%hold on
h1=num2str(age);

xlim([t_start t_max+t_future])


hold on
plot([t_max,time_future],[M_wLC(age+1,end),M_wLC_Forecast(age+1,:)],'r--')
%plot(time_future,M_wLC_Forecast_UB(age+1,:),'r:')
%plot(time_future,M_wLC_Forecast_LB(age+1,:),'r:')
hold on
plot([t_max,time_future],[M_D(age+1,end),M_D_Forecast(age+1,:)],'k--')
%plot(time_future,M_D_Forecast_UB(age+1,:),'k:')
%plot(time_future,M_D_Forecast_LB(age+1,:),'k:')
hold on
plot([t_max,time_future],[M_ML(age+1,end),M_ML_Forecast(age+1,:)],'b--')
%plot(time_future,M_ML_Forecast_UB(age+1,:),'b:')
%plot(time_future,M_ML_Forecast_LB(age+1,:),'b:')
hold on
%plot([t_max,time_future],[M_W(age+1,end),M_W_Forecast(age+1,:)],'c--')
%plot(time_future,M_W_Forecast_UB(age+1,:),'c:')
%plot(time_future,M_W_Forecast_LB(age+1,:),'c:')
%legend('historical','wLC','D','ML','W')
%legend('historical','wLC','D','ML')

age=65; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose

tempM=M_full(age+1,:);
index42=tempM~=42;

plot(time_full(index42),tempM(index42),'g')
h=num2str(age);
title(ax1,['m(x,t) over time for age ' h1 ' (steeper) and ' h])
xlim([t_start t_max+t_future])

hold on
plot([t_max,time_future],[M_wLC(age+1,end),M_wLC_Forecast(age+1,:)],'r--')
hold on
plot([t_max,time_future],[M_D(age+1,end),M_D_Forecast(age+1,:)],'k--')
hold on
plot([t_max,time_future],[M_ML(age+1,end),M_ML_Forecast(age+1,:)],'b--')
hold on
%plot([t_max,time_future],[M_W(age+1,end),M_W_Forecast(age+1,:)],'c--')


ax1 = subplot(4,1,2);
%plots development ver time for a specific age
age=40; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose

tempM=M_full(age+1,:);
index42=tempM~=42;

plot(time_full(index42),tempM(index42),'g')

age1=age;
xlim([t_start t_max+t_future])

hold on
plot([t_max,time_future],[M_wLC(age+1,end),M_wLC_Forecast(age+1,:)],'r--')
hold on
plot([t_max,time_future],[M_D(age+1,end),M_D_Forecast(age+1,:)],'k--')
hold on
plot([t_max,time_future],[M_ML(age+1,end),M_ML_Forecast(age+1,:)],'b--')
hold on
%plot([t_max,time_future],[M_W(age+1,end),M_W_Forecast(age+1,:)],'c--')

age=20; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose

tempM=M_full(age+1,:);
index42=tempM~=42;

plot(time_full(index42),tempM(index42),'g')
h=num2str(age);
h1=num2str(age1);

title(ax1,['m(x,t) over time for age ' h ' (lower) and age ' h1])
xlim([t_start t_max+t_future])

hold on
plot([t_max,time_future],[M_wLC(age+1,end),M_wLC_Forecast(age+1,:)],'r--')
hold on
plot([t_max,time_future],[M_D(age+1,end),M_D_Forecast(age+1,:)],'k--')
hold on
plot([t_max,time_future],[M_ML(age+1,end),M_ML_Forecast(age+1,:)],'b--')
hold on
%plot([t_max,time_future],[M_W(age+1,end),M_W_Forecast(age+1,:)],'c--')

ax1 = subplot(4,1,3);
%plots development over time for a specific age
age=99; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose

tempM=M_full(age+1,:);
index42=tempM~=42;

plot(time_full(index42),tempM(index42),'g')
h=num2str(age);
title(ax1,['m(x,t) over time for age ' h])
xlim([t_start t_max+t_future])

hold on
plot([t_max,time_future],[M_wLC(age+1,end),M_wLC_Forecast(age+1,:)],'r--')
hold on
plot([t_max,time_future],[M_D(age+1,end),M_D_Forecast(age+1,:)],'k--')
hold on
plot([t_max,time_future],[M_ML(age+1,end),M_ML_Forecast(age+1,:)],'b--')
hold on
%plot([t_max,time_future],[M_W(age+1,end),M_W_Forecast(age+1,:)],'c--')
ax1 = subplot(4,1,4);
%plots development over time for a specific age
age=108; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose

tempM=M_full(age+1,:);
index42=tempM~=42;

plot(time_full(index42),tempM(index42),'g')
h=num2str(age);
title(ax1,['m(x,t) over time for age ' h])
xlim([t_start t_max+t_future])

hold on
plot([t_max,time_future],[M_wLC(age+1,end),M_wLC_Forecast(age+1,:)],'r--')
hold on
plot([t_max,time_future],[M_D(age+1,end),M_D_Forecast(age+1,:)],'k--')
hold on
%plot([t_max,time_future],[M_ML(age+1,end),M_ML_Forecast(age+1,:)],'b--')
hold on
%plot([t_max,time_future],[M_W(age+1,end),M_W_Forecast(age+1,:)],'c--')
set(forecast1,'Units','Inches');
pos = get(forecast1,'Position');
set(forecast1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(forecast1,'forecast1','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forecast2=figure;
%plots development over age for a specific year
year=19; %beetween 1:t_future) only forecasted years %%%%choose
if t+year<=t_full
    templnM=lnM_full(ages+1,t+year);
    indexln42=templnM~=log(42);
    plot(ages(indexln42),templnM(indexln42),'g')
    hold on
end
plot(ages,lnM_wLC_Forecast(:,year),'r')
hold on
plot(ages,lnM_D_Forecast(:,year),'k')
hold on
plot(ages,lnM_ML_Forecast(:,year),'b')
%hold on
%plot(ages,lnM_W_Forecast(:,year),'c')
h=num2str(year+t_max);
h1= num2str( ages(1));
h2= num2str(ages(end));
title(['ln[m(x,t)] for ages ' h1 ' to ' h2 ' in year ' h])

%hold on
%plot(ages,lnM_wLC_Forecast_LB(:,year),'r:')
%plot(ages,lnM_wLC_Forecast_UB(:,year),'r:')
%hold on
%plot(ages,lnM_D_Forecast_UB(:,year),'k:')
%plot(ages,lnM_D_Forecast_LB(:,year),'k:')
%hold on
%plot(ages,lnM_ML_Forecast_UB(:,year),'b:')
%plot(ages,lnM_ML_Forecast_LB(:,year),'b:')
%hold on
%plot(ages,lnM_W_Forecast_UB(:,year),'c:')
%plot(ages,lnM_W_Forecast_LB(:,year),'c:')

%year=t_future; %beetween 1:t_future) only forecasted years %%%%choose
%if t+year<=t_full
%    templnM=lnM_full(ages+1,t+year);
%    indexln42=templnM~=log(42);
%    plot(ages(indexln42),templnM(indexln42),'g')
%    hold on
%end
%plot(ages,lnM_wLC_Forecast(:,year),'r')
%hold on
%plot(ages,lnM_D_Forecast(:,year),'k')
%hold on
%plot(ages,lnM_ML_Forecast(:,year),'b')
%hold on
%plot(ages,lnM_W_Forecast(:,year),'c')
%}
%year=1; %beetween 1:t_future) only forecasted years %%%%choose
%if t+year<=t_full
%    templnM=lnM_full(ages+1,t+year);
%    indexln42=templnM~=log(42);
%    plot(ages(indexln42),templnM(indexln42),'g')
%    hold on
%end
%plot(ages,lnM_wLC_Forecast(:,year),'r')
%hold on
%plot(ages,lnM_D_Forecast(:,year),'k')
%hold on
%plot(ages,lnM_ML_Forecast(:,year),'b')
%hold on
%plot(ages,lnM_W_Forecast(:,year),'c')
legend({'historical','wLC','D','ML'},'Location','northwest')
set(forecast2,'Units','Inches');
pos = get(forecast2,'Position');
set(forecast2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(forecast2,'forecast2','-dpdf','-r0')

%%%%%%%%%%
errforecast=figure;
subplot(2,2,1)
tempM=M_full(:,t+1:end);
index42=tempM==42;
Merr= tempM- M_ML_Forecast(:,1:t_full-t);
Merr(index42)=0;
xscale = [1997 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
%xscale = [1996 2015];

yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
%xscale = [1996 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
%xscale = [1996 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
suptitle('ML forecast error colormap')

set(errforecast,'Units','Inches');
pos = get(errforecast,'Position');
set(errforecast,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errforecast,'errforecastML','-dpdf','-r0')

errforecast=figure;
subplot(2,2,1)
tempM=M_full(:,t+1:end);
index42=tempM==42;
Merr= tempM- M_wLC_Forecast(:,1:t_full-t);
Merr(index42)=0;
xscale = [1997 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
%xscale = [1996 2015];

yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
%xscale = [1996 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
%xscale = [1996 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
suptitle('1-0 weights forecast error colormap')

set(errforecast,'Units','Inches');
pos = get(errforecast,'Position');
set(errforecast,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errforecast,'errforecastwLC','-dpdf','-r0')

errforecast=figure;
subplot(2,2,1)
tempM=M_full(:,t+1:end);
index42=tempM==42;
Merr= tempM- M_D2_Forecast(:,1:t_full-t);
Merr(index42)=0;
xscale = [1997 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
%xscale = [1996 2015];

yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
%xscale = [1996 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
%xscale = [1996 2015];
yscale = [0 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
suptitle('Death weights forecast error colormap ')

set(errforecast,'Units','Inches');
pos = get(errforecast,'Position');
set(errforecast,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errforecast,'errforecastD','-dpdf','-r0')

b=figure;
plot(linspace(1,x,x),b_wLC,'r');
hold on;
plot(linspace(1,x,x),b_D,'k');
hold on;
plot(linspace(1,x,x),b_ML,'b');
title('b_x')
legend({'wLC','D','ML'},'Location','northwest')

set(b,'Units','Inches');
pos = get(b,'Position');
set(b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(b,'bhistory','-dpdf','-r0')



clearvars h h1 h2 n age ages ax1 ax2 i j time year time_future lnM_SVD_Forecast lnM_SVD_Forecast_UB
clearvars lnM_SVD_Forecast_LB lnM_ML_Forecast lnM_ML_Forecast_UB lnM_ML_Forecast_LB lnM_LC_Forecast
clearvars lnM_LC_Forecast_UB lnM_LC_Forecast_LB lnM_full lnM_SVD lnM_ML lnM_LC 
clearvars index42 tempM time_full indexln42 templnM