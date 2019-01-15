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

%below function were from gdp model, set g_t to 0 so we can use the
%functions for the simple case fitting


%no weights but second step adjusted
tic
[a_wLC,b_wLC,k_wLC,steps_wLC,fval_wLC]=Fit_wLC(lnM);
time_wLC=toc;
tic
[a_newton,b_newton,k_newton,fval_newton,exitflag_newton2]=FitLC_newton(lnM,zeros(x,1),ones(x,1),ones(t,1));
time_newton=toc;
[a_simplex,b_simplex,k_simplex,fval_simplex,exitflag_simplex2]=FitLC_simplex(lnM,zeros(x,1),ones(x,1),ones(t,1));

[a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);
[a_newton2,b_newton2,k_newton2]=LeeCarterFit(D,E,a_newton,b_newton,k_newton);
[a_simplex2,b_simplex2,k_simplex2]=LeeCarterFit(D,E,a_simplex,b_simplex,k_simplex);
 
%weights given by D
tic
[a_DwLC,b_DwLC,k_DwLC,steps_D,fval_D]=Fit_wLC(lnM,D);
timeDwLC=toc;
tic
[a_Dnewton,b_Dnewton,k_Dnewton,fval_Dnewton,exitflag_Dnewton2]=FitLC_newton(lnM,zeros(x,1),ones(x,1),ones(t,1),D);
[a_Dnewton,b_Dnewton,k_Dnewton,fval_Dnewton,exitflag_Dnewton2]=FitLC_simplex(lnM,a_Dnewton,b_Dnewton,k_Dnewton,D);
[a_Dnewton,b_Dnewton,k_Dnewton,fval_Dnewton,exitflag_Dnewton2]=FitLC_newton(lnM,a_Dnewton,b_Dnewton,k_Dnewton,D);
time_Dnewton=toc
 tic
 [a_Dsimplex,b_Dsimplex,k_Dsimplex,fval_Dsimplex,exitflag_Dsimplex]=FitLC_simplex(lnM,zeros(x,1),ones(x,1),ones(t,1),D);
 time_Dsimplex=toc;

%custom weights
%w_x,t=t+x_max-x 
W=zeros(x,t);
for i=0:x_max
    for j=1:t
    W(i+1,j)=t_min+j-1+x_max-i;
    end
end
 B=lnM==log(42);
    W(B)=0;
    
 [a_W,b_W,k_W,steps_W,fval_W]=Fit_wLC(lnM,W);
 [a_Wnewton,b_Wnewton,k_Wnewton,fval_Wnewton,exitflag_Wnewton]=FitLC_newton(lnM,zeros(x,1),ones(x,1),ones(t,1),W);
 [a_Wsimplex,b_Wsimplex,k_Wsimplex,fval_Wsimplex,exitflag_Wsimplex]=FitLC_simplex(lnM,zeros(x,1),ones(x,1),ones(t,1),W);

 [a_W2,b_W2,k_W2]=LeeCarterFit(D,E,a_W,b_W,k_W);