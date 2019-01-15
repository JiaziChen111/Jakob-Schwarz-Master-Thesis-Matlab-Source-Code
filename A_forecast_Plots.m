%calc alpha/2 and 1-alpha/2 quantil path for A_1 and A_2 sample path
    %alpha=0.1; set in Master
    alpha=0.1;
    [A_1_qlower,A_1_qupper,A_2_qlower,A_2_qupper]=QuantilPath(A_future,alpha);
    [A_1__unc_qlower,A_1_unc_qupper,A_2_unc_qlower,A_2_unc_qupper]=QuantilPath(A_future_uncertain,alpha);
    

data_time=linspace(t_min_C,t_max_C,t_C);
forecast_time=linspace(t_max_C+1,t_max_C+t_future_C,t_future_C);
forecast_time=[forecast_time(1)-1,forecast_time];
view=1:t_C; %%%choose start year to view past fit
Atforecast2=figure('name','A(t) Development and Forecast');
ax1 = subplot(2,1,1);
%plots development over time for A1
plot(data_time(view),A_1(view),'g');
hold on
plot(forecast_time,[A_1(t_C);A_future_mean(:,1)],'r');
hold on
plot(forecast_time,[A_1(t_C);A_1_qupper],'--');
hold on
plot(forecast_time,[A_1(t_C);A_1_unc_qupper],':');
hold on
plot(forecast_time,[A_1(t_C);A_1_qlower],'--');
hold on
plot(forecast_time,[A_1(t_C);A_1__unc_qlower],':');
title(ax1,'A1(t) over time with forecast')

ax2 = subplot(2,1,2);
%plots development over time for A1
plot(data_time(view),A_2(view),'g');
hold on
plot(forecast_time,[A_2(t_C);A_future_mean(:,2)],'r');
hold on
plot(forecast_time,[A_2(t_C);A_2_qupper],'--');
hold on
plot(forecast_time,[A_2(t_C);A_2_unc_qupper],':');
hold on
plot(forecast_time,[A_2(t_C);A_2_qlower],'--');
hold on
plot(forecast_time,[A_2(t_C);A_2_unc_qlower],':');
title(ax2,'A2(t) over time with forecast')

legend({'historical fit','mean forecast','90% quantil', 'with parameter uncertainty'},'Location','Best')


set(Atforecast2,'Units','Inches');
pos = get(Atforecast2,'Position');
set(Atforecast2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(Atforecast2,'Atforecast2M','-dpdf','-r0')

%clearvars data_time forecast_time view ax1 ax2 A_1_unc_qupper A_2_unc_qlower
%clearvars A_1_qlower A_1_qupper A_2_qlower A_2_qupper A_1__unc_qlower A_2_unc_qupper
