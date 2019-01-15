clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=100; %new maximum age group 100+, set to prevent 0 and NaN in M
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
        
%%%Fitting Plots, no forecasting%%%
    %M_Plots;
    %q_Plots;
    %M_CombinedPlots;
    
%%%Forecasting%%%
    
    %Box_Jenkins; %plots a_x b_x k_t A(t) dA(t) dk(t) and autocorrealtions and partial autoco. for k and dk
    
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
    t_future=19;
    n_paths=100;
    
    %run Arima
    
    A_x=[a_wLC,a_wLC2,a_W,a_W2,a_D,a_D2,a_ML,a_ML2];
     B_x=[b_wLC,b_wLC2,b_W,b_W2,b_D,b_D2,b_ML,b_ML2];
      K_t=[k_wLC,k_wLC2,k_W,k_W2,k_D,k_D2,k_ML,k_ML2];
      
      ErrM=zeros(9*8,6);
      Err2Mat=zeros(9,8);
      MeanErrMat=zeros(9,8);
      
    for i=1:8
        a_x=A_x(:,i);
        b_x=B_x(:,i);
        k_t=K_t(:,i);
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima010,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+1,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(1,i)=Err2;
    MeanErrMat(1,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima011,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+2,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(2,i)=Err2;
    MeanErrMat(2,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima110,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+3,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(3,i)=Err2;
    MeanErrMat(3,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima111,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+4,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(4,i)=Err2;
    MeanErrMat(4,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima012,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+5,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(5,i)=Err2;
    MeanErrMat(5,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima112,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+6,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(6,i)=Err2;
    MeanErrMat(6,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima210,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+7,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(7,i)=Err2;
    MeanErrMat(7,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima211,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+8,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(8,i)=Err2;
    MeanErrMat(8,i)=MeanErr;
    
    [Est_LC,EstCov_LC,res_LC,Mean_LC,UB_LC,LB_LC, Ysim_LC ,simU_LC, simL_LC, YsimParUnc_LC, U_ParUnc_LC, L_ParUnc_LC] = ArimaForecast(t_future,n_paths,Arima212,k_t);  
   [ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ]  = errorfkt( M_full(x_min+1:x_max+1,t+1:end),D_full(x_min+1:x_max+1,t+1:end),E_full(x_min+1:x_max+1,t+1:end),a_x,b_x,Mean_LC );
    ErrM(8*(i-1)+9,:)=[MaxErr,MinErr,MeanErr,Err2,L2,FVU];
    Err2Mat(9,i)=Err2;
    MeanErrMat(9,i)=MeanErr;
    
    end 
    
  fprintf([repmat('%.5f\t', 1, size(Err2Mat, 2)) '\n'], Err2Mat')
   fprintf([repmat('%.5f\t', 1, size(Err2Mat, 2)) '\n'], MeanErrMat')
    
    
    
    
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
    
    %M_ForecastPlots;