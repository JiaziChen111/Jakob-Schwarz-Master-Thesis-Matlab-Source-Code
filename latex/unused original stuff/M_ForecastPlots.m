n=size(M);
x=n(1);
t=n(2);

M_LC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_LC(i,j)=exp(a_LC(i)+b_LC(i)*k_LC(j));
        end
    end

M_ML=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML(i,j)=exp(a_ML(i)+b_ML(i)*k_ML(j));
        end
    end
   
M_SVD=zeros(x,t);
    for i=1:x
        for j=1:t
            M_SVD(i,j)=exp(a_SVD(i)+b_SVD(i)*k_SVD(j));
        end
    end
    
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
    
lnM_SVD = log(M_SVD);
lnM_ML = log(M_ML);
lnM_LC = log(M_LC);
lnM_wLC = log(M_wLC);
lnM_wLC2 = log(M_wLC2);

M_LC_Forecast=zeros(x,t_future);
M_LC_Forecast_UB=zeros(x,t_future);
M_LC_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_LC_Forecast(i,j)=exp(a_LC(i)+b_LC(i)*Mean_LC(j));
             M_LC_Forecast_UB(i,j)=exp(a_LC(i)+b_LC(i)*UB_LC(j));
              M_LC_Forecast_LB(i,j)=exp(a_LC(i)+b_LC(i)*LB_LC(j));
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
   
M_SVD_Forecast=zeros(x,t_future);
M_SVD_Forecast_UB=zeros(x,t_future);
M_SVD_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
           M_SVD_Forecast(i,j)=exp(a_SVD(i)+b_SVD(i)*Mean_SVD(j));
            M_SVD_Forecast_UB(i,j)=exp(a_SVD(i)+b_SVD(i)*UB_SVD(j));
              M_SVD_Forecast_LB(i,j)=exp(a_SVD(i)+b_SVD(i)*LB_SVD(j));
        end
    end
    
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

lnM_SVD_Forecast=log(M_SVD_Forecast);
lnM_SVD_Forecast_UB=log(M_SVD_Forecast);
lnM_SVD_Forecast_LB=log(M_SVD_Forecast);

lnM_ML_Forecast=log(M_ML_Forecast);
lnM_ML_Forecast_UB=log(M_ML_Forecast_UB);
lnM_ML_Forecast_LB=log(M_ML_Forecast_LB);

lnM_LC_Forecast=log(M_LC_Forecast);
lnM_LC_Forecast_UB=log(M_LC_Forecast_UB);
lnM_LC_Forecast_LB=log(M_LC_Forecast_LB);

lnM_wLC_Forecast=log(M_wLC_Forecast);
lnM_wLC_Forecast_UB=log(M_wLC_Forecast_UB);
lnM_wLC_Forecast_LB=log(M_wLC_Forecast_LB);

lnM_wLC2_Forecast=log(M_wLC2_Forecast);
lnM_wLC2_Forecast_UB=log(M_wLC2_Forecast_UB);
lnM_wLC2_Forecast_LB=log(M_wLC2_Forecast_LB);

lnM_full=log(M_full);    
    
    
figure('name','m(x,t) Fits and Forecast of Lee-Carter Modell');
ax1 = subplot(2,1,1);
%plots development over time for a specific age
age=40; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

tempM=M_full(age,:);
index42=tempM~=42;
time_full=linspace(t_start,t_end,t_full);

plot(time_full(index42),tempM(index42),'g')
hold on
plot(time,M_LC(age,:),'r')
hold on
plot(time,M_SVD(age,:),'k')
hold on
plot(time,M_ML(age,:),'b')
hold on
plot(time,M_wLC(age,:),'c')
hold on
plot(time,M_wLC2(age,:),'m')
h=num2str(ages(age));
title(ax1,['m(x,t) over time for age ' h])
legend('historical','LC','SVD','ML','wLC','wLC2')

time_future=linspace(t_max+1,t_max+t_future,t_future);
hold on
plot([t_max,time_future],[M_LC(age,end),M_LC_Forecast(age,:)],'r--')
plot(time_future,M_LC_Forecast_UB(age,:),'r:')
plot(time_future,M_LC_Forecast_LB(age,:),'r:')
hold on
plot([t_max,time_future],[M_SVD(age,end),M_SVD_Forecast(age,:)],'k--')
plot(time_future,M_SVD_Forecast_UB(age,:),'k:')
plot(time_future,M_SVD_Forecast_LB(age,:),'k:')
hold on
plot([t_max,time_future],[M_ML(age,end),M_ML_Forecast(age,:)],'b--')
plot(time_future,M_ML_Forecast_UB(age,:),'b:')
plot(time_future,M_ML_Forecast_LB(age,:),'b:')
hold on
plot([t_max,time_future],[M_wLC(age,end),M_wLC_Forecast(age,:)],'c--')
plot(time_future,M_wLC_Forecast_UB(age,:),'c:')
plot(time_future,M_wLC_Forecast_LB(age,:),'c:')
hold on
plot([t_max,time_future],[M_wLC2(age,end),M_wLC2_Forecast(age,:)],'m--')
plot(time_future,M_wLC2_Forecast_UB(age,:),'m:')
plot(time_future,M_wLC2_Forecast_LB(age,:),'m:')

ax2 = subplot(2,1,2);
%plots development over age for a specific year
year=19; %beetween 1:t_future) only forecasted years %%%%choose
if t+year<=t_full
    templnM=lnM_full(ages+1,t+year);
    indexln42=templnM~=log(42);
    plot(ages(indexln42),templnM(indexln42),'g')
    hold on
end
plot(ages,lnM_LC_Forecast(:,year),'r')
hold on
plot(ages,lnM_SVD_Forecast(:,year),'k')
hold on
plot(ages,lnM_ML_Forecast(:,year),'b')
hold on
plot(ages,lnM_wLC_Forecast(:,year),'c')
hold on
plot(ages,lnM_wLC2_Forecast(:,year),'m')
h=num2str(year+t_max);
h1= num2str( ages(1));
h2= num2str(ages(end));
title(ax2,['ln[m(x,t)] for ages ' h1 ' to ' h2 ' in year ' h])

hold on
plot(ages,lnM_LC_Forecast_LB(:,year),'r:')
plot(ages,lnM_LC_Forecast_UB(:,year),'r:')
hold on
plot(ages,lnM_SVD_Forecast_UB(:,year),'k:')
plot(ages,lnM_SVD_Forecast_LB(:,year),'k:')
hold on
plot(ages,lnM_ML_Forecast_UB(:,year),'b:')
plot(ages,lnM_ML_Forecast_LB(:,year),'b:')
hold on
plot(ages,lnM_wLC_Forecast_UB(:,year),'c:')
plot(ages,lnM_wLC_Forecast_LB(:,year),'c:')

hold on
plot(ages,lnM_wLC2_Forecast_UB(:,year),'m:')
plot(ages,lnM_wLC2_Forecast_LB(:,year),'m:')


clearvars h h1 h2 n age ages ax1 ax2 i j time year time_future lnM_SVD_Forecast lnM_SVD_Forecast_UB
clearvars lnM_SVD_Forecast_LB lnM_ML_Forecast lnM_ML_Forecast_UB lnM_ML_Forecast_LB lnM_LC_Forecast
clearvars lnM_LC_Forecast_UB lnM_LC_Forecast_LB lnM_full lnM_SVD lnM_ML lnM_LC 
clearvars index42 tempM time_full indexln42 templnM