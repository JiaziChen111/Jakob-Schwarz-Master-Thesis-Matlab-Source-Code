clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=110; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=1996; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

%load data with above parameters
BuildCentralDeathMatrix;

    %%%Cairns model
   x_max_C=110; %for a new maximum age group, i .e. 100 instead of 110
   x_min_C=70; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
        x_max_C_70=110; %for a new maximum age group, i .e. 100 instead of 110
        x_min_C_70= 70; %for a new minimum included age group, i .e. 60 instead of 0
        x_C_70=x_C;
        [A_1_70, A_2_70,mu_70,C_70]=q_RegressionFit(q,x_min_C_70,x_max_C_70);

    

q_C_70=zeros(x_C,t_C);
    for i=1:x_C
        for j=1:t_C
            h=A_1_70(j)+A_2_70(j)*((x_min_C+i-1)+(j-1));
            q_C_70(i,j)=exp(h)/(1+exp(h));
        end
    end

    
%%%Cairns model
   x_max_C=110; %for a new maximum age group, i .e. 100 instead of 110
   x_min_C= 40; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1, A_2,mu,C]=q_RegressionFit(q,x_min_C,x_max_C);

    q_C=zeros(x_C,t_C);
    for i=1:x_C
        for j=1:t_C
            h=A_1(j)+A_2(j)*((x_min_C+i-1)+(j-1));
            q_C(i,j)=exp(h)/(1+exp(h));
        end
    end
    

%no weights but second step adjusted
%use implied M by q
lnM=log(-log(1-q_full(:,t_min-1955:t_min-1955+t-1)));
[a_wLC,b_wLC,k_wLC,steps_wLC,fval_wLC]=Fit_wLC(lnM);
  
%%%%%%%
M_wLC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC(i,j)=exp(a_wLC(i)+b_wLC(i)*k_wLC(j));
        end
    end
 
    
q_wLC=1-exp(-M_wLC);

 
    arima211=arima(2,1,1);
    
    %Parameters for Arima Forecast and Plot
    t_future=24;
    n_paths=10;
   [Est_wLC,~,res_wLC,Mean_wLC,UB_wLC,LB_wLC, Ysim_wLC ,simU_wLC, simL_wLC, YsimParUnc_wLC, U_ParUnc_wLC, L_ParUnc_wLC] = ArimaForecast(t_future,n_paths,arima211,k_wLC);  
    M_wLC_211=zeros(x,t_future);
    %M_wLC_Forecast_UB=zeros(x,t_future);
    %M_wLC_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_wLC_211(i,j)=exp(a_wLC(i)+b_wLC(i)*Mean_wLC(j));
   %          M_wLC_Forecast_UB(i,j)=exp(a_wLC(i)+b_wLC(i)*UB_wLC(j));
   %           M_wLC_Forecast_LB(i,j)=exp(a_wLC(i)+b_wLC(i)*LB_wLC(j));
        end
    end
    q_wLC_211=1-exp(-M_wLC_211);
   % q_wLC_Forecast_UB=1-exp(-M_wLC_Forecast_UB);
   % q_wLC_Forecast_LB=1-exp(-M_wLC_Forecast_LB);
arima011=arima(0,1,1);
     [Est_wLC,~,res_wLC,Mean_wLC,UB_wLC,LB_wLC, Ysim_wLC ,simU_wLC, simL_wLC, YsimParUnc_wLC, U_ParUnc_wLC, L_ParUnc_wLC] = ArimaForecast(t_future,n_paths,arima011,k_wLC);  
    M_wLC_011=zeros(x,t_future);
    %M_wLC_Forecast_UB=zeros(x,t_future);
    %M_wLC_Forecast_LB=zeros(x,t_future);
    for i=1:x
        for j=1:t_future
            M_wLC_011(i,j)=exp(a_wLC(i)+b_wLC(i)*Mean_wLC(j));
   %          M_wLC_Forecast_UB(i,j)=exp(a_wLC(i)+b_wLC(i)*UB_wLC(j));
   %           M_wLC_Forecast_LB(i,j)=exp(a_wLC(i)+b_wLC(i)*LB_wLC(j));
        end
    end
    q_wLC_011=1-exp(-M_wLC_011);
   % q_wLC_Forecast_UB=1-exp(-M_wLC_Forecast_UB);
   % q_wLC_Forecast_LB=1-exp(-M_wLC_Forecast_LB)
   
   
   
   %%%Forecasting At%%%
    
    %Box_Jenkins; %plots a_x b_x k_t A(t) dA(t) dk(t) and autocorrealtions and partial autoco. for k and dk

    %Forecasting A(t) by Cairns model including parameter uncertainty
    t_future_C=24; %number of periods we want to forecast into the future 
    n_paths_C=10; %number of paths
    %Creates expeted future path of A
    A_future_mean=A_MeanPath(A_1(t),A_2(t),mu,t_future_C);
    %creates a cell array size n_paths of simulated paths length t_future
    A_future=A_SamplePaths(A_1(t),A_2(t),mu,C,t_future_C,n_paths_C);
    
    %parameter relizations from the posteriori dist given back as cell array
    %[mu_uncertain,C_uncertain]=Normal_Inverse_Wishart(A_1,A_2,n_paths_C);
    [mu_hat,V_hat,df]=PosterioriParameter(A_1,A_2);
    %creates a cell array size n_paths, each patch uses a diffeent mu and C
    A_future_uncertain=A_UncertainSamplePaths(A_1(t),A_2(t),mu_hat,V_hat,df,t_future_C,n_paths_C);
    clearvars mu_hat V_hat df
    
    %A_forecast_Plots;
    %q_ForecastPlots;
    
    %calc alpha/2 and 1-alpha/2 quantil for q_C from sample path
    alpha=0.1;
    q_C_future_UB=zeros(x_C,t_future_C);
    q_C_future_LB=zeros(x_C,t_future_C);
    H=zeros(n_paths_C,1);
    for i=1:x_C
        for j=1:t_future_C
            for r=1:n_paths_C
                h=A_future{r}(j,1)+A_future{r}(j,2)*((x_min_C+i-1)+(j-1+t_C));
                H(r,1)=exp(h)/(1+exp(h));
            end
       q_C_future_LB(i,j) = prctile(H',100*(alpha/2),2);
       q_C_future_UB(i,j) = prctile(H',100*(1-alpha/2),2);    
        end
    end

    lnq_C_future_LB=log(q_C_future_LB);
    lnq_C_future_UB=log(q_C_future_UB);
    
    Dt=[diff(A_1),diff(A_2)];
    Dfuture=cell(6,1);
    for i=0:5
        Md=varm(2,i);
        EstMd = estimate(Md,Dt);
        isstable(EstMd)
        Dfuture{i+1} = forecast(EstMd,t_future_C,Dt);
    end
    Afuture=cell(6,1);
    for j=1:6
        Afuture{j}=zeros(t_future_C,2);
        Afuture{j}(1,:)=Dfuture{j}(i,:)+[A_1(end),A_2(end)];
        for i=1:t_future_C-1
            Afuture{j}(i+1,:)=Dfuture{j}(i+1,:)+Afuture{j}(i,:);
        end
    end
    q_C_future=cell(6,1);
    for k=1:6
        q_C_future{k}=zeros(x_C,t_future_C);
        for i=1:x_C
            for j=1:t_future_C    
                    h=Afuture{k}(j,1)+Afuture{k}(j,2)*((x_min_C+i-1)+(j+t_C-1));
                    q_C_future{k}(i,j)=exp(h)/(1+exp(h));
            end
        end
    end
     Dt_70=[diff(A_1_70),diff(A_2_70)];
       Dfuture_70=cell(6,1);
    for i=0:5
        Md=varm(2,i);
        EstMd = estimate(Md,Dt_70);
        isstable(EstMd)
        Dfuture_70{i+1} = forecast(EstMd,t_future_C,Dt_70);
    end
    Afuture_70=cell(6,1);
    for j=1:6
        Afuture_70{j}=zeros(t_future_C,2);
        Afuture_70{j}(1,:)=Dfuture_70{j}(i,:)+[A_1_70(end),A_2_70(end)];
        for i=1:t_future_C-1
            Afuture_70{j}(i+1,:)=Dfuture_70{j}(i+1,:)+Afuture_70{j}(i,:);
        end
    end
    q_C_future_70=cell(6,1);
   for k=1:6
        q_C_future_70{k}=zeros(x_C_70,t_future_C);
        for i=1:x_C_70
            for j=1:t_future_C    
                    h=Afuture_70{k}(j,1)+Afuture_70{k}(j,2)*((x_min_C_70+i-1)+(j+t_C-1));
                    q_C_future_70{k}(i,j)=exp(h)/(1+exp(h));
            end
        end
    end
 q_true=q_full(x_min_C+1:x_max_C+1,(t_max_C-t_start_C+2):t_full_C);
%%%%%%%%%%%%%%
k=1;
a1=figure;
subplot(4,1,1)
age=40;

h=num2str(age);
plot(linspace(1956,2015,60),q_full(age+1,:),'g')
hold on;
plot(linspace(1956,1996,41),q_wLC(age+1,:),'c')
hold on
plot(linspace(1956,1996,41),q_C(age-x_min_C+1,:),'b')
hold on;
plot(linspace(1996,2020,25),[q_wLC(age+1,end),q_wLC_011(age+1,:)],'c--')
hold on;
plot(linspace(1996,2020,25),[q_C(age-x_min_C+1,end),q_C_future{k}(age-x_min_C+1,:)],'b--')

%legend('historical','wLC011','Cairns','Cairns_{70}')
title(['q(x,t) over time for age ' h])

subplot(4,1,2)
age=70;
h=num2str(age);
plot(linspace(1956,2015,60),q_full(age+1,:),'g')
hold on;
plot(linspace(1956,1996,41),q_wLC(age+1,:),'c')
hold on
plot(linspace(1956,1996,41),q_C(age-x_min_C+1,:),'b')
hold on;
plot(linspace(1956,1996,41),q_C_70(age-x_min_C_70+1,:),'r')
hold on;
plot(linspace(1996,2020,25),[q_wLC(age+1,end),q_wLC_011(age+1,:)],'c--')
hold on;
plot(linspace(1996,2020,25),[q_C(age-x_min_C+1,end),q_C_future{k}(age-x_min_C+1,:)],'b--')
hold on;
plot(linspace(1996,2020,25),[q_C_70(age-x_min_C_70+1,end),q_C_future_70{k}(age-x_min_C_70+1,:)],'r--')

%legend('historical','wLC011','Cairns','Cairns_{70}')
title(['q(x,t) over time for age ' h])

subplot(4,1,3)
age=95;
h=num2str(age);
plot(linspace(1956,2015,60),q_full(age+1,:),'g')
hold on;
plot(linspace(1956,1996,41),q_wLC(age+1,:),'c')
hold on
plot(linspace(1956,1996,41),q_C(age-x_min_C+1,:),'b')
hold on;
plot(linspace(1956,1996,41),q_C_70(age-x_min_C_70+1,:),'r')
hold on;
plot(linspace(1996,2020,25),[q_wLC(age+1,end),q_wLC_011(age+1,:)],'c--')
hold on;
plot(linspace(1996,2020,25),[q_C(age-x_min_C+1,end),q_C_future{k}(age-x_min_C+1,:)],'b--')
hold on;
plot(linspace(1996,2020,25),[q_C_70(age-x_min_C_70+1,end),q_C_future_70{k}(age-x_min_C_70+1,:)],'r--')

%legend('historical','wLC011','Cairns','Cairns_{70}')
title(['q(x,t) over time for age ' h])

subplot(4,1,4)
age=110;
h=num2str(age);
plot(linspace(1956,2015,60),q_full(age+1,:),'g')
hold on;
plot(linspace(1956,1996,41),q_wLC(age+1,:),'c')
hold on
plot(linspace(1956,1996,41),q_C(age-x_min_C+1,:),'b')
hold on;
plot(linspace(1956,1996,41),q_C_70(age-x_min_C_70+1,:),'r')
hold on;
plot(linspace(1996,2020,25),[q_wLC(age+1,end),q_wLC_011(age+1,:)],'c--')
hold on;
plot(linspace(1996,2020,25),[q_C(age-x_min_C+1,end),q_C_future{k}(age-x_min_C+1,:)],'b--')
hold on;
plot(linspace(1996,2020,25),[q_C_70(age-x_min_C_70+1,end),q_C_future_70{k}(age-x_min_C_70+1,:)],'r--')

%legend('historical','wLC011','Cairns','Cairns_{70}')
title(['q(x,t) over time for age ' h])

set(a1,'Units','Inches');
pos = get(a1,'Position');
set(a1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(a1,'qfuture','-dpdf','-r0')

a2=figure;
h=num2str(age);
plot(linspace(40,110,71),log(q_true(:,end)),'g')
hold on;
plot(linspace(40,110,71),log(q_wLC_011(41:end,19)),'c')
hold on;
plot(linspace(40,110,71),log(q_C_future{1}(:,19)),'b')
hold on;
plot(linspace(70,110,41),log(q_C_future_70{1}(:,19)),'r')


legend('historical','wLC011','Cairns','Cairns_{70}','location','best')
title('ln[q(x,t)] forecast over age in 2015')

set(a2,'Units','Inches');
pos = get(a2,'Position');
set(a2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(a2,'qfuture2','-dpdf','-r0')


%%%%%%%%%%
errq=figure;
subplot(2,2,1)
Merr= q_true- q_C_future{1}(:,1:19) ;

xscale = [1997 2015];
yscale = [40 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.3 0.3]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.05]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,4)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.001 0.001]);
xlabel('year t')
ylabel('age x')
title('')

suptitle('Fitting error colormap: Cairns')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errqfuture','-dpdf','-r0')

%%%%%%%%%%
errq=figure;
subplot(2,2,1)
Merr= q_true- q_wLC(41,1:19) ;

xscale = [1997 2015];
yscale = [40 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.3 0.3]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.05]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,4)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.001 0.001]);
xlabel('year t')
ylabel('age x')
title('')

suptitle('Fitting error colormap: LC')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errqfuturewLC','-dpdf','-r0')

%%%%%%%%%%
errq=figure;
subplot(2,1,1)
Merr= q_true(31:end,:)- q_C_future_70{1}(:,1:19) ;

xscale = [1997 2015];
yscale = [70 110];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.3 0.3]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,1,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')
title('')

suptitle('Fitting error colormap: Cairns70')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errqfuture70','-dpdf','-r0')