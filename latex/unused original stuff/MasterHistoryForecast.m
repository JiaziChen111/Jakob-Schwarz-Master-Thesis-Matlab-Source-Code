clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=104; %new maximum age group 100+, set to prevent 0 and NaN in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for fit
t_max=1996; %last year of data included for fit
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=0; %do not replace 0 and -1 with 42 in M, instead here we set x_max=104

%load data with above parameters
BuildCentralDeathMatrix; 


%%%Lee Carter Model Fitting%%%
       
    %calc SVD least sqarres solution for fit
    [a_SVD,b_SVD,k_SVD]=SVDFit(lnM);
    [a_wLC,b_wLC,k_wLC]=Fit_wLC(lnM,D);
    %k_t adjustment to replicate observed total deaths with estimated m(x,t)
    [a_LC,b_LC,k_LC]=LeeCarterFit(D,E,a_SVD,b_SVD,k_SVD);
    [a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);
    %note that sum k != 0 after reestimation
    
    %following code plots the error remaining in the reestimated k_LC
    %{
    tt=linspace(t_min,t_max,t);
    err=zeros(1,t);
    for j=1:t
        for i=1:x
            err(j)=err(j)+E(i,j)*exp((a_LC(i)+b_LC(i)*k_LC(j)));
        end
        err(j)=abs(sum (D(:,j))-err(j)); %sum of all deaths year t-SumE);
    end
    plot(tt,err)
    title('Sum(x) D(x,t)=Sum(x) N(x,t)*e^(a_x+k_t*b_x)')
    legend('absolut error')
    max(err) %0.0018 for males in west germany data 1958-2015 0-100+
    sum(err) %.0142
    sqrt(sum(err.^2)) % 0.0029
    %}
    
%%%Poisson modeling Fitting%%%
   
    [a_ML,b_ML,k_ML]=PoissonFit(D,E);    
    
%%%Cairns Model%%%
    x_max_C=109; %for a new maximum age group, i .e. 100 instead of 110
    x_min_C= 40; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=1956; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=1996; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1, A_2,mu,C]=q_RegressionFit(q,x_min_C,x_max_C);

%%%Fitting Plots, no forecasting%%%
    %M_Plots;
    %q_Plots;
    %M_CombinedPlots;
    
%%%Forecasting%%%
    
    %Box_Jenkins; %plots a_x b_x k_t A(t) dA(t) dk(t) and autocorrealtions and partial autoco. for k and dk

    %Forecasting A(t) by Cairns model including parameter uncertainty
    t_future_C=19; %number of periods we want to forecast into the future 1997-2015
    n_paths_C=1000; %number of paths
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
    q_ForecastPlots;
    
    Arima010=arima(0,1,0);
    Arima011=arima(0,1,1);
    
    %Parameters for Arima Forecast and Plot
    t_future=19;
    n_paths=100;
    
    %run Arima
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima011,k_LC);
    [Est_SVD,EstCov_SVD,res_SVD,Mean_SVD,UB_SVD,LB_SVD, Ysim_SVD ,simU_SVD, simL_SVD, YsimParUnc_SVD, U_ParUnc_SVD, L_ParUnc_SVD] = ArimaForecast(t_future,n_paths,Arima011,k_SVD);
    [Est_ML,EstCov_ML,res_ML,Mean_ML,UB_ML,LB_ML, Ysim_ML ,simU_ML, simL_ML, YsimParUnc_ML, U_ParUnc_ML, L_ParUnc_ML] = ArimaForecast(t_future,n_paths,Arima011,k_ML);
    [Est_wLC,EstCov_wLC,res_wLC,Mean_wLC,UB_wLC,LB_wLC, Ysim_wLC ,simU_wLC, simL_wLC, YsimParUnc_wLC, U_ParUnc_wLC, L_ParUnc_wLC] = ArimaForecast(t_future,n_paths,Arima011,k_wLC);
    [Est_wLC2,EstCov_wLC2,res_wLC2,Mean_wLC2,UB_wLC2,LB_wLC2, Ysim_wLC2 ,simU_wLC2, simL_wLC2, YsimParUnc_wLC2, U_ParUnc_wLC2, L_ParUnc_wLC2] = ArimaForecast(t_future,n_paths,Arima011,k_wLC2);
   
    %{
    %choose what to plot
    yF=Mean_LC;
    res = res_LC;
    Est=Est_LC;
    kt=k_LC;
    EstCov=EstCov_LC;
    fit= 'Lee Carter';
    model= 'Arima(0,1,0)';    
    UB=UB_LC;
    LB=LB_LC;
    simU=simU_LC;
    U_ParUnc=U_ParUnc_LC;
    L_ParUnc=L_ParUnc_LC;
    simL=simL_LC;
   
    Arima_Plots;
    %}
    
    M_ForecastPlots;