clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=105; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1970; %first year of data included for LC like models
t_max=2005; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

%load data with above parameters
BuildCentralDeathMatrix; 

BuildGDPvector;
BuildDrugData;




  W=ones(x,t);
  H=lnM==log(42);
  W(H)=0;

%cut out of 1951-2017 data
%heatdays_cut=heatdays(t_min-1950:t_min-1950+t-1);
%cut out of 1951-2018 
%heatnights_cut=heatnights(t_min-1950:t_min-1950+t-1);
%cut inflation of out 1950-2017 data
%inflation_cut=inflation(t_min-1949:t_min-1949+t-1);
%cut out of alkohol 1970-2005
alk_cut=alk(t_min-1969:t_min-1969+t-1);
%cut out of 1950-2017 (2006-2017 estimated) data alk_2
%alk_2_cut=alk_2(t_min-1949:t_min-1949+t-1);
%cut out of 1979-2015 data
%raum_f_young_cut=rauch_f_young(t_min-1978:t_min-1978+t-1);
%cut out of 1950-2003 data
%rauch_f_cut=rauch_f(t_min-1949:t_min-1949+t-1);
%cut out of 1964-2017 data
zig_cut=zig(t_min-1963:t_min-1963+t-1);

g=[g_t';zig_cut';alk_cut'];
n=size(g);
n=n(1);

tic
[a_ex_newton,b_ex_newton,k_ex_newton,y_ex_newton,fval_ex_newton,exitflag_ex_newton]=FitLC_Explanation_newton(lnM,g,zeros(x,1),ones(x,1),ones(t,1),zeros(x,n),W,0);
timeLC_newton=toc;

%tic   
%[a_ML_newton,b_ML_newton,k_ML_newton,y_ML_newton,fval_ML_newton,exitflag_ML_newton]=PoissonFit_explanation_newton(D, E,g, zeros(x,1),ones(x,1),zeros(t,1),zeros(x,n),0);
%timeML_newton=toc;

tic
[a_ML,b_ML,k_ML,y_ML,fval_ML]=PoissonFit_explanation(D, E,g, zeros(x,1),ones(x,1),zeros(t,1),zeros(x,n),0);
timeML=toc;

%[a_x,b_x,k_t,fval,exitflag,steps]=Fit_wLC(lnM,D);
%tic
%[a_normalen,b_normalen,k_normalen,y_normalen,fval_normalen,exitflag_normalen,steps_normalen]=FitLC_GDP_normalen(lnM,g_t,zeros(x,1),ones(x,1),ones(t,1),ones(x,1),W);
%timeLC_normalen=toc;
%tic
%[a_ex_newton,b_ex_newton,k_ex_newton,y_ex_newton,fval_ex_newton,exitflag_ex_newton]=FitLC_Explanation_newton(lnM,g_t',zeros(x,1),ones(x,1),ones(t,1),zeros(x,1),W,0);
%timeLC_newton=toc;
%tic
%[a_gdp_newton,b_gdp_newton,k_gdp_newton,y_gdp_newton,fval_gdp_newton,exitflag_gdp_newton]=FitLC_GDP_newton(lnM,g_t,zeros(x,1),ones(x,1),ones(t,1),zeros(x,1),W);
%timeLC_newton_gdp=toc;

hhh=figure;
subplot(2,2,1);
plot(linspace(1,x,x),a_ex_newton(1:x),'r');
hold on;
%plot(linspace(1,x,x),a_ML_newton(1:x),'b');
%hold on;
plot(linspace(1,x,x),a_ML(1:x),'b');
title('a_x');

subplot(2,2,2);
plot(linspace(1,x,x),b_ex_newton,'r');
hold on;
%plot(linspace(1,x,x),b_ML_newton,'b');
%hold on;
plot(linspace(1,x,x),b_ML,'b');
title('b_x');

subplot(2,2,3:4);
plot(linspace(1,t,t),k_ex_newton,'r');
%hold on;
%plot(linspace(1,t,t),k_ML_newton,'b');
hold on;
plot(linspace(1,t,t),k_ML,'b');

title('k_t');
%legend('wLC_{newton}','ML_{newton}','ML_{1dim}');
legend('wLC_{newton}','ML_{1dim}');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(hhh,'inflationabk','-dpdf','-r0')

hhhhh=figure;
subplot(2,1,1)
plot(linspace(1,t,t),g(1,:))
title('log GDP/capita');

subplot(2,1,2);
plot(linspace(1,x,x),y_ex_newton(:,1),'r');
hold on;
%plot(linspace(1,x,x),y_ML_newton(:,1),'b');
%hold on;
plot(linspace(1,x,x),y_ML(:,1),'b');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

hhhh=figure;
subplot(2,2,1)
plot(linspace(1,t,t),g(3,:))
title('alkohol per capita');

subplot(2,2,2);
plot(linspace(1,x,x),y_ex_newton(:,3),'r');
hold on;
%plot(linspace(1,x,x),y_ML_newton(:,1),'b');
%hold on;
plot(linspace(1,x,x),y_ML(:,3),'b');

subplot(2,2,3)
plot(linspace(1,t,t),g(2,:))
title('sold zigarettes');

subplot(2,2,4);
plot(linspace(1,x,x),y_ex_newton(:,2),'r');
hold on;
%plot(linspace(1,x,x),y_ML_newton(:,2),'b');
%hold on;
plot(linspace(1,x,x),y_ML(:,2),'b');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(hhhh,'inflationy','-dpdf','-r0')