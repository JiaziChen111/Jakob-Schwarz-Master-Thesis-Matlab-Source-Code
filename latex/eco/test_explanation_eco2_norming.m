clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=110; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

%load data with above parameters
BuildCentralDeathMatrix; 

BuildGDPvector;
%BuildDrugData;

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
%alk_cut=alk(t_min-1969:t_min-1969+t-1);
%cut out of 1950-2017 (2006-2017 estimated) data alk_2
%alk_2_cut=alk_2(t_min-1949:t_min-1949+t-1);
%cut out of 1979-2015 data
%raum_f_young_cut=rauch_f_young(t_min-1978:t_min-1978+t-1);
%cut out of 1950-2003 data
%rauch_f_cut=rauch_f(t_min-1949:t_min-1949+t-1);
%cut out of 1964-2017 data
%zig_cut=zig(t_min-1963:t_min-1963+t-1);

g=[g_t'];
n=size(g);
n=n(1);
norming=1;

a_start=zeros(x,1);
b_start=ones(x,1);
k_start=ones(t,1);
y_start=zeros(x,n);

a_startML=zeros(x,1);
b_startML=ones(x,1);
k_startML=zeros(t,1);
y_startML=zeros(x,n);



tic
[a_2step,b_2step,k_2step,y_2step,fval_2step,exitflag_2step,steps_2steps]=FitLC_Explanation_2step(lnM,g,norming,W);
time2step=toc;

tic
[a_ex_newton,b_ex_newton,k_ex_newton,y_ex_newton,fval_ex_newton,exitflag_ex_newton]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,W,norming);
timeLC_newton=toc;
%tic
%[a_ex_newton2,b_ex_newton2,k_ex_newton2,y_ex_newton2,~,~]=FitLC_Explanation_simplex(lnM,g,a_ex_newton,b_ex_newton,k_ex_newton,y_ex_newton,W,norming);
%[a_ex_newton2,b_ex_newton2,k_ex_newton2,y_ex_newton2,fval_ex_newton2,exitflag_ex_newton2]=FitLC_Explanation_newton(lnM,g,a_ex_newton2,b_ex_newton2,k_ex_newton2,y_ex_newton2,W,norming);
%timeLC_newton2=toc;

tic   
[a_ML_newton,b_ML_newton,k_ML_newton,y_ML_newton,fval_ML_newton,exitflag_ML_newton]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,norming);
timeML_newton=toc;

FUN = @(z) sum(sum(W.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)'-reshape(z(2*x+t+1:end),x,n)*g).^2   ));
        y0=reshape(y_ML_newton,x*n,1);
        x0 = [a_ML_newton;b_ML_newton;k_ML_newton;y0];
        fval2_ML_newton=FUN(x0);

tic
[a_ML,b_ML,k_ML,y_ML,fval_ML]=PoissonFit_explanation(D, E,g, a_startML,b_startML,k_startML,y_startML,norming);
timeML=toc;

        y0=reshape(y_ML,x*n,1);
        x0 = [a_ML;b_ML;k_ML;y0];
        fval2_ML=FUN(x0);


%tic
%[a_normalen,b_normalen,k_normalen,y_normalen,fval_normalen,exitflag_normalen,steps_normalen]=FitLC_GDP_normalen(lnM,g_t,a_start,b_start,k_start,y_start,W);
%timeLC_normalen=toc;

a_start=a_2step;
b_start=b_2step;
k_start=k_2step;
y_start=y_2step;

a_startML=a_2step;
b_startML=b_2step;
k_startML=k_2step;
y_startML=y_2step;

tic   
[a_ML_newton2,b_ML_newton2,k_ML_newton2,y_ML_newton2,fval_ML_newton2,exitflag_ML_newton2]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,norming);
timeML_newton2=toc;

        y0=reshape(y_ML_newton2,x*n,1);
        x0 = [a_ML_newton2;b_ML_newton2;k_ML_newton2;y0];
        fval2_ML_newton2=FUN(x0);

tic
[a_ML2,b_ML2,k_ML2,y_ML2,fval_ML2]=PoissonFit_explanation(D, E,g, a_startML,b_startML,k_startML,y_startML,norming);
timeML2=toc;

        y0=reshape(y_ML2,x*n,1);
        x0 = [a_ML2;b_ML2;k_ML2;y0];
        fval2_ML2=FUN(x0);

  tic
[a_ex_newton2,b_ex_newton2,k_ex_newton2,y_ex_newton2,fval_ex_newton2,exitflag_ex_newton2]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,W,norming);
timeLC_newton2=toc;

x_cut=101   ;   
hhh=figure;
subplot(2,2,1);
plot(linspace(1,x_cut ,x_cut ),a_ex_newton2(1:x_cut ),'k');
hold on;
%plot(linspace(1,x_cut ,x_cut ),a_ex_newton(1:x_cut ),'r');
hold on;
plot(linspace(1,x_cut ,x_cut ),a_ML_newton2(1:x_cut ),'c');
hold on;
plot(linspace(1,x_cut ,x_cut ),a_ML(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut ,x_cut ),a_2step(1:x_cut ),'y');
hold on;
%plot(linspace(1,x_cut ,x_cut ),a_normalen(1:x_cut ),'g');

title('a_x');

subplot(2,2,2);
plot(linspace(1,x_cut ,x_cut ),b_ex_newton2(1:x_cut ),'k');
hold on;
%plot(linspace(1,x_cut ,x_cut ),b_ex_newton(1:x_cut ),'r');
hold on;
plot(linspace(1,x_cut ,x_cut ),b_ML_newton2(1:x_cut ),'c');
hold on;
plot(linspace(1,x_cut ,x_cut ),b_ML(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut ,x_cut ),b_2step(1:x_cut ),'y');
hold on;
%plot(linspace(1,x_cut ,x_cut ),b_normalen(1:x_cut ),'g');


title('b_x');


subplot(2,2,3:4);
plot(linspace(1,t,t),k_ex_newton2(1:t),'k');
hold on;
%plot(linspace(1,t,t),k_ex_newton(1:t),'r');
hold on;
plot(linspace(1,t,t),k_ML_newton2(1:t),'c');
hold on;
plot(linspace(1,t,t),k_ML(1:t),'b');
hold on;
plot(linspace(1,t,t),k_2step(1:t),'y');
hold on;
%plot(linspace(1,t,t),k_normalen(1:t),'g');


title('k_t');
legend({'wLC_{newton3step}','ML_{newton}','ML_{1dim}','3step'},'location','best');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh,'gdpabknorming','-dpdf','-r0')

hhhhh=figure;

plot(linspace(1,x_cut ,x_cut ),y_ex_newton2(1:x_cut ),'k');
hold on;
%plot(linspace(1,x_cut ,x_cut ),y_ex_newton(1:x_cut ),'r');
hold on;
plot(linspace(1,x_cut ,x_cut ),y_ML_newton2(1:x_cut ),'c');
hold on;
plot(linspace(1,x_cut ,x_cut ),y_ML(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut ,x_cut ),y_2step(1:x_cut ),'y');
hold on;
%plot(linspace(1,x_cut ,x_cut ),y_normalen(1:x_cut ),'g');



title('y_x');
legend({'wLC_{newton3step}','ML_{newton}','ML_{1dim}','3step'},'location','best');
set(hhhhh,'Units','Inches');
pos = get(hhhhh,'Position');
set(hhhhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhhh,'gdpynorming','-dpdf','-r0')

