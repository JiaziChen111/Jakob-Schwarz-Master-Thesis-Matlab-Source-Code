clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=100; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

%load data with above parameters
BuildCentralDeathMatrix; 


%%%Lee Carter Model Fitting%%%
    
    %{
    %calc SVD least sqarres solution for fit
    [a_LC,b_LC,k_SVD]=SVDFit(lnM);    
    [a_wLC,b_wLC,k_wLC]=Fit_wLC(lnM);
    %for xmin=0 xmax=103 and all 60 years we get
    norm(a_LC-a_wLC)+norm(b_LC-b_wLC)+norm(k_SVD-k_wLC)
    %=6.7547e-08
    %very small difference. wLC has the advatage of that we can avoid the
    zero cell problem with weights, see wilmorth 1993
    %}
   
    [a_SVD,b_SVD,k_SVD]=Fit_wLC(lnM);
    [a_wLC,b_wLC,k_wLC]=Fit_wLC(lnM,D);
    %norm(a_LC-a_wLC)+norm(b_LC-b_wLC)+norm(k_SVD-k_wLC) 
    %for age 0-109 all years, very small difference in a and b, high for k
    
    %k_t adjustment to replicate observed total deaths with estimated m(x,t)
     [a_LC,b_LC,k_LC]=LeeCarterFit(D,E,a_SVD,b_SVD,k_SVD);
    [a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);
    %norm(k_LC-k_wLC2) smaller difference but still high
    %note that sum k != 0 after reestimation
    
    %following code plots the error remaining in the reestimated k_LC
    %{
    tt=linspace(t_min,t_max,t);
    err=zeros(1,t);
    for j=1:t
        for i=1:x
            err(j)=err(j)+E(i,j)*exp((a_wLC(i)+b_wLC(i)*k_wLC2(j)));
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
     %the LC solutions input is optional as a start value for the iteration
    [a_ML,b_ML,k_ML]=PoissonFit(D,E,a_LC,b_LC,k_LC);    
    
%%%Cairns Model%%%
    x_max_C=110; %for a new maximum age group, i .e. 100 instead of 110
    x_min_C= 40; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=1956; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=2015; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1, A_2,mu,C]=q_RegressionFit(q,x_min_C,x_max_C);

%%%Fitting Plots, no forecasting%%%
    M_Plots;
    q_Plots;
    M_CombinedPlots;
    
%%%Forecasting%%%
    
    Box_Jenkins; %plots a_x b_x k_t A(t) dA(t) dk(t) and autocorrealtions and partial autoco. for k and dk

    %Forecasting A(t) by Cairns model including parameter uncertainty
    t_future_C=110; %number of periods we want to forecast into the future 
    n_paths_C=5000; %number of paths
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
    
    alpha=0.1;
    A_forecast_Plots;
    q_ForecastPlots;
    
    Arima010=arima(0,1,0);
    Arima011=arima(0,1,1);
    
    %Parameters for Arima Forecast and Plot
    t_future=95;
    n_paths=100;
    
    %run Arima
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima011,k_LC);
    [Est_SVD,EstCov_SVD,res_SVD,Mean_SVD,UB_SVD,LB_SVD, Ysim_SVD ,simU_SVD, simL_SVD, YsimParUnc_SVD, U_ParUnc_SVD, L_ParUnc_SVD] = ArimaForecast(t_future,n_paths,Arima011,k_SVD);
    [Est_ML,EstCov_ML,res_ML,Mean_ML,UB_ML,LB_ML, Ysim_ML ,simU_ML, simL_ML, YsimParUnc_ML, U_ParUnc_ML, L_ParUnc_ML] = ArimaForecast(t_future,n_paths,Arima011,k_ML);
    [Est_wLC,EstCov_wLC,res_wLC,Mean_wLC,UB_wLC,LB_wLC, Ysim_wLC ,simU_wLC, simL_wLC, YsimParUnc_wLC, U_ParUnc_wLC, L_ParUnc_wLC] = ArimaForecast(t_future,n_paths,Arima011,k_wLC);
    [Est_wLC2,EstCov_wLC2,res_wLC2,Mean_wLC2,UB_wLC2,LB_wLC2, Ysim_wLC2 ,simU_wLC2, simL_wLC2, YsimParUnc_wLC2, U_ParUnc_wLC2, L_ParUnc_wLC2] = ArimaForecast(t_future,n_paths,Arima011,k_wLC2);
   
    Arima_Plots;
    %%%%
    M_ForecastPlots;
    
    %%%Remaining Life Expectantcy Calculation with [alpha/2,1-alpha/2] CI
    %[e_C] = lifeexp_q([1 26 46], [5 15 25 35 45], [q,q_C_future],x_min_C, x_max_C);
    [e_C,e_C_UB,e_C_LB,qsim] = lifeexp_q([1 26 46], [61], [q,q_C_future],x_min_C, x_max_C, A_future, alpha);
    
    %create a Cell array, each cell is a realized M matrix fit+forecast
    Msim_LC=getMsim(M_LC,a_LC,b_LC,Ysim_LC); %used for CI below for e_xt
    [e_LC, e_LC_UB, e_LC_LB] = lifeexp_m([1 41 66 86], [5 15 25 35 45], [M_LC,M_LC_Forecast],x_min, x_max, Msim_LC,alpha);
    [e_ML] = lifeexp_m([1 41 66 86], [5 15 25 35 45], [M_ML,M_ML_Forecast],x_min, x_max);
    [e_SVD] = lifeexp_m([1 41 66 86], [5 15 25 35 45], [M_SVD,M_SVD_Forecast],x_min, x_max);
    [e_wLC] = lifeexp_m([1 41 66 86], [5 15 25 35 45], [M_wLC,M_wLC_Forecast],x_min, x_max);
    [e_wLC2] = lifeexp_m([1 41 66 86], [5 15 25 35 45], [M_wLC2,M_wLC2_Forecast],x_min, x_max);
    
    
%%%GDP-Model%%%
    BuildGDPvector;
    %first we fit the model: ln m(x,t)=y_0x+y_1x*g_t+e_x,t
    [y_0x,y_1x]=PureGDPFit(lnM,g_t); 
    
    %fit the extended GDP Lee Carter model: ln m(x,t)=a_x+b_x*k_t+y_x*g_t+e
    [a_gdp_newton,b_gdp_newton,k_gdp_newton,y_gdp_newton,fval_newton,exitflag_newton]=FitLC_GDP_newton(lnM,g_t,a_SVD,b_SVD,k_SVD,zeros(x,1));
    [a_gdp_newton2,b_gdp_newton2,k_gdp_newton2,y_gdp_newton2,fval_newton2,exitflag_newton2]=FitLC_GDP_newton(lnM,g_t,ones(x,1),ones(x,1)/x,zeros(t,1),zeros(x,1));
    %same solution for both starting values, lowest sqarred error
    
    [a_gdp_simplex,b_gdp_simplex,k_gdp_simplex,y_gdp_simplex,fval_simplex,exitflag_simplex]=FitLC_GDP_simplex(lnM,g_t,a_SVD,b_SVD,k_SVD,zeros(x,1));
    %seems to very slowly converge to newton solution, better than normalen
    %[a_gdp_simplex2,b_gdp_simplex2,k_gdp_simplex2,y_gdp_simplex2,fval_simplex2,exitflag_simplex2]=FitLC_GDP_simplex(lnM,g_t,ones(x,1),ones(x,1)/x,ones(t,1),ones(x,1));
    %doesnt work, too far away
   
    [a_gdp_normalen,b_gdp_normalen,k_gdp_normalen,y_gdp_normalen,fval_normalen,exitflag_normalen]=FitLC_GDP_normalen(lnM,g_t,a_gdp_newton,b_gdp_newton,k_gdp_newton,y_gdp_newton);
    [a_gdp_normalen2,b_gdp_normalen2,k_gdp_normalen2,y_gdp_normalen2,fval_normalen2,exitflag_normalen2]=FitLC_GDP_normalen(lnM,g_t,ones(x,1),ones(x,1)/x,ones(t,1),ones(x,1));
    %same solution for both starting values. close to newton solution 
    %but worse, doesnt change any further no matter the tolerace settings
    
    %solutions can be improved by alterating methods
    %[a_gdp_simplex,b_gdp_simplex,k_gdp_simplex,y_gdp_simplex,fval_simplex,exitflag_simplex]=FitLC_GDP_simplex(lnM,g_t,a_gdp_newton,b_gdp_newton,k_gdp_newton,y_gdp_newton,D);
    %[a_gdp_newton,b_gdp_newton,k_gdp_newton,y_gdp_newton,fval_newton,exitflag_newton]=FitLC_GDP_newton(lnM,g_t,a_gdp_simplex,b_gdp_simplex,k_gdp_simplex,y_gdp_simplex,D);
    %[a_gdp_simplex,b_gdp_simplex,k_gdp_simplex,y_gdp_simplex,fval_simplex,exitflag_simplex]=FitLC_GDP_simplex(lnM,g_t,a_gdp_newton,b_gdp_newton,k_gdp_newton,y_gdp_newton,D);
    
    
    
    plot_gdp;