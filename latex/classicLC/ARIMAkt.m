clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=104; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

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

dk_wLC=diff(k_wLC);
ddk_wLC=diff(diff(k_wLC));

dk=figure('name','Identifying d');
subplot(3,2,1);
t=length(k_wLC);
plot(linspace(1,t,t),k_wLC,'r');
title('k_{wLC}');

subplot(3,2,2) 
autocorr(k_wLC)
ylabel(' ')
title('Sample Autocorr.Function for k_{wLC}')

subplot(3,2,3) 
h=length(dk_wLC);
plot(linspace(1,h,h),dk_wLC,'r');
title('dk_{wLC}');
subplot(3,2,4)
autocorr(dk_wLC)
ylabel(' ')
title('Sample Autocorr. Function for dk_{wLC}')

subplot(3,2,5) 
h=length(ddk_wLC);
plot(linspace(1,h,h),ddk_wLC,'r');
title('ddk_{wLC}')
subplot(3,2,6)
autocorr(ddk_wLC)
ylabel(' ')
title('Sample Autocorr. Function for ddk_{wLC}')

set(dk,'Units','Inches');
pos = get(dk,'Position');
set(dk,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(dk,'dk','-dpdf','-r0')
    
pq=figure('name','ACF and PACF for dk_t');
dk_wLC=diff(k_wLC);
dk_wLC2=diff(k_wLC2);
dk_ML=diff(k_ML);
dk_D=diff(k_D);
dk_W=diff(k_W);
subplot(3,2,[1 2]);
plot(linspace(1,t-1,t-1),dk_wLC,'r');
hold on;
plot(linspace(1,t-1,t-1),dk_ML,'b');
title('dk_t=k_t-k_{t-1}');

subplot(3,2,3) 
autocorr(dk_wLC)
ylabel(' ')
title('Sample ACF for dk_{wLC}')
subplot(3,2,4)
parcorr(dk_wLC)
ylabel(' ')
title('Sample PACF for dk_{wLC}')

subplot(3,2,5) 
autocorr(dk_ML)
ylabel(' ')
title('Sample ACF for dk_{ML}')
subplot(3,2,6)
parcorr(dk_ML)
ylabel(' ')
title('Sample PACF for dk_{ML}')

set(pq,'Units','Inches');
pos = get(pq,'Position');
set(pq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(pq,'pq2','-dpdf','-r0')

 
%%calc AICBIC table
AIC=zeros(9,6);
BIC=zeros(9,6);
Kt=[k_wLC,k_wLC2,k_W,k_W2,k_ML,k_D];
h=size(Kt);
for i=1:h(2)
kt=Kt(:,i);

Model=arima(0,1,0);
[Est,EstCov,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(1,i)=aic;
BIC(1,i)=bic;

Model=arima(0,1,1);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(2,i)=aic;
BIC(2,i)=bic;

Model=arima(1,1,0);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(3,i)=aic;
BIC(3,i)=bic;

Model=arima(1,1,1);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(4,i)=aic;
BIC(4,i)=bic;

Model=arima(0,1,2);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(5,i)=aic;
BIC(5,i)=bic;

Model=arima(1,1,2);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(6,i)=aic;
BIC(6,i)=bic;

Model=arima(2,1,0);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(7,i)=aic;
BIC(7,i)=bic;

Model=arima(2,1,1);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(8,i)=aic;
BIC(8,i)=bic;

Model=arima(2,1,2);
[~,~,logL,~]=estimate(Model,kt);
parameter=2+length(Model.MA)+length(Model.AR);
[aic,bic] = aicbic(logL, parameter, length(kt));
AIC(9,i)=aic;
BIC(9,i)=bic;

end
%{


%%%Forecasting%%%
    Arima010=arima(0,1,0);
    Arima011=arima(0,1,1);
    Arima110=arima(1,1,0);
    Arima111=arima(1,1,1);
    Arima012=arima(0,1,2);
    Arima112=arima(1,1,2);
    Arima210=arima(2,1,0);
    Arima211=arima(2,1,1);
    Arima212=arima(2,1,2);
    
    
    %Parameters for Arima Forecast and Plot
    t_future=50;
    n_paths=100;
    
    %run Arima
    tt=figure;
    subplot(2,1,1)
    t=length(k_wLC);
    plot(linspace(t_min,t_max,t),k_wLC,'g');
    
    title('k_t (females)')
    [~,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC,aic_LC,bic_LC] = ArimaForecast(t_future,n_paths,Arima010,k_wLC);
    hold on;
    plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);Mean_LC],'b-')
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);UB_LC],'b--');
    hasbehavior(h, 'legend', false); 
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);LB_LC],'b--');
    hasbehavior(h, 'legend', false); 
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC,aic_LC,bic_LC] = ArimaForecast(t_future,n_paths,Arima011,k_wLC);
     hold on;
    plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);Mean_LC],'r-')
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);UB_LC],'r--');
    hasbehavior(h, 'legend', false); 
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);LB_LC],'r--');
    hasbehavior(h, 'legend', false); 
   
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, ~, ~, U_ParUnc_LC, L_ParUnc_LC,aic_LC,bic_LC] = ArimaForecast(t_future,n_paths,Arima211,k_wLC);
     hold on;
    plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);Mean_LC],'k-')
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);UB_LC],'k--');
    hasbehavior(h, 'legend', false); 
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);LB_LC],'k--');
    hasbehavior(h, 'legend', false); 
    
    legend({'k_{wLC}','arima010','arima011','arima211'},'Location','southwest');
    
    subplot(2,1,2)
    %%%%%%%
    DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =4; % chose sex 3 4 5 for female male total
x_max=104; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

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
    
 h=plot(linspace(t_min,t_max,t),k_wLC,'g');
    hasbehavior(h, 'legend', false); 
    title('k_t (males)')
    [~,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC,aic_LC,bic_LC] = ArimaForecast(t_future,n_paths,Arima010,k_wLC);
    hold on;
    plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);Mean_LC],'b-')
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);UB_LC],'b--');
    hasbehavior(h, 'legend', false); 
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);LB_LC],'b--');
    hasbehavior(h, 'legend', false); 
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC,aic_LC,bic_LC] = ArimaForecast(t_future,n_paths,Arima011,k_wLC);
     hold on;
    plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);Mean_LC],'r-')
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);UB_LC],'r--');
    hasbehavior(h, 'legend', false); 
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);LB_LC],'r--');
    hasbehavior(h, 'legend', false); 
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, ~, ~, U_ParUnc_LC, L_ParUnc_LC,aic_LC,bic_LC] = ArimaForecast(t_future,n_paths,Arima211,k_wLC);
     hold on;
    plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);Mean_LC],'k-')
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);UB_LC],'k--');
    hasbehavior(h, 'legend', false); 
    hold on;
    h=plot(linspace(t_max,t_max+t_future,t_future+1),[k_wLC(end);LB_LC],'k--');
    hasbehavior(h, 'legend', false); 
    
    
    
    
    
    %%%%%
    set(tt,'Units','Inches');
pos = get(tt,'Position');
set(tt,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(tt,'arimas','-dpdf','-r0')

    %[Est_SVD,EstCov_SVD,res_SVD,Mean_SVD,UB_SVD,LB_SVD, Ysim_SVD ,simU_SVD, simL_SVD, YsimParUnc_SVD, U_ParUnc_SVD, L_ParUnc_SVD] = ArimaForecast(t_future,n_paths,Arima011,k_SVD);
    %[Est_ML,EstCov_ML,res_ML,Mean_ML,UB_ML,LB_ML, Ysim_ML ,simU_ML, simL_ML, YsimParUnc_ML, U_ParUnc_ML, L_ParUnc_ML] = ArimaForecast(t_future,n_paths,Arima011,k_ML);
    %[Est_wLC,EstCov_wLC,res_wLC,Mean_wLC,UB_wLC,LB_wLC, Ysim_wLC ,simU_wLC, simL_wLC, YsimParUnc_wLC, U_ParUnc_wLC, L_ParUnc_wLC] = ArimaForecast(t_future,n_paths,Arima011,k_wLC);
    %[Est_wLC2,EstCov_wLC2,res_wLC2,Mean_wLC2,UB_wLC2,LB_wLC2, Ysim_wLC2 ,simU_wLC2, simL_wLC2, YsimParUnc_wLC2, U_ParUnc_wLC2, L_ParUnc_wLC2] = ArimaForecast(t_future,n_paths,Arima011,k_wLC2);
   %}

