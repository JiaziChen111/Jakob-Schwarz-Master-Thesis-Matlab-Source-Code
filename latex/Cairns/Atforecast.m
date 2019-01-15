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
        x_min_C_70= 40; %for a new minimum included age group, i .e. 60 instead of 0
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
    t_future=19;
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
    t_future_C=19; %number of periods we want to forecast into the future 
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
    
    A_forecast_Plots;
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
        Dfuture{i+1} = forecast(EstMd,t_future_C,Dt);
    end
    Afuture=cell(6,1);
    for j=1:6
        Afuture{j}=zeros(19,2);
        Afuture{j}(1,:)=Dfuture{j}(i,:)+[A_1(end),A_2(end)];
        for i=1:t_future_C-1
            Afuture{j}(i+1,:)=Dfuture{j}(i+1,:)+Afuture{j}(i,:);
        end
    end
    q_C_future=cell(6,1);
    for k=1:6
        q_C_future{k}=zeros(x_C,t_future_C);
        for i=1:x_C
            for j=1:t_C+t_future_C
                if j<=t_C
                    continue
                else            
                    h=Afuture{k}(j-t_C,1)+Afuture{k}(j-t_C,2)*((x_min_C+i-1)+(j-1));
                    q_C_future{k}(i,j-t_C)=exp(h)/(1+exp(h));
                end
            end
        end
    end
 q_true=q_full(x_min_C_70+1:x_max_C_70+1,(t_max_C-t_start_C+2):t_full_C);
 
 ErrM=zeros(8,2);
 for i=1:6
    [ ~,~,MeanErr,Err2,~,~]=errorfkt_q( q_true,q_C_future{i}(1:end,:) );
    ErrM(i,1)=Err2;
 end
 [ ~,~,MeanErr,Err2,~,~]=errorfkt_q( q_true, q_wLC_011(x_min_C_70+1:x_max_C_70+1,:) );
 ErrM(7,1)=Err2;
 [ ~,~,MeanErr,Err2,~,~]=errorfkt_q( q_true, q_wLC_211(x_min_C_70+1:x_max_C_70+1,:) );
 ErrM(8,1)=Err2;
  
 
 
 
fprintf([repmat('%.5f\t', 1, size(ErrM, 2)) '\n'], ErrM')
