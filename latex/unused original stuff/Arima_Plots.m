 yF=Mean_LC;
    res = res_LC;
    Est=Est_LC;
    kt=k_LC;
    EstCov=EstCov_LC;
    fit= 'LC';
    model= 'Arima(0,1,0)';    
    UB=UB_LC;
    LB=LB_LC;
    simU=simU_LC;
    U_ParUnc=U_ParUnc_LC;
    L_ParUnc=L_ParUnc_LC;
    simL=simL_LC;

figure('Name',['normalized residuals of k_' fit])
subplot(2,2,1)
plot(res./sqrt(Est.Variance))
title('Standardized Residuals')
subplot(2,2,2)
qqplot(res)
subplot(2,2,3)
autocorr(res)
subplot(2,2,4)
parcorr(res)

figure
h4 = plot(kt,'Color',[.75,.75,.75]);
hold on
h5 = plot(t:t+t_future,[kt(end);yF],'r','LineWidth',2);

h6 = plot(t+1:t+t_future,UB,'Color',[.35,.85,.85],'LineWidth',2);
plot(t+1:t+t_future,LB,'Color',[.35,.85,.85],'LineWidth',2);

h7 = plot(t+1:t+t_future,simU,'--','LineWidth',1.5);
plot(t+1:t+t_future,simL,'--','LineWidth',1.5);

h8 = plot(t+1:t+t_future,U_ParUnc,':','LineWidth',1.5);
plot(t+1:t+t_future,L_ParUnc,':','LineWidth',1.5);

legend([h4,h5 ,h6,h7,h8],'historical','Forecast',...
       'Forecast Interval 95%','Forecast MC 95%','MC including ParUnc','Location','Northwest')
title(['k_{' fit '} Forecast with ' model])
hold off

clearvars kt Est UB LB simU simL yMSE res fit model h4 h5 h6 h7 temp mu_temp temp_Ysim EstCov U_ParUnc
