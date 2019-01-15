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

%no weights but second step adjusted

[a_wLC,b_wLC,k_wLC,steps_wLC,fval_wLC]=Fit_wLC(lnM);
[a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);

%weights given by D
[a_DwLC,b_DwLC,k_DwLC,steps_D,fval_D]=Fit_wLC(lnM,D);
[a_DwLC2,b_DwLC2,k_DwLC2]=LeeCarterFit(D,E,a_DwLC,b_DwLC,k_DwLC);

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
 [a_W2,b_W2,k_W2]=LeeCarterFit(D,E,a_W,b_W,k_W);
 
  [a_ML,b_ML,k_ML,fval_ML]=PoissonFit(D, E);
  [a_ML2,b_ML2,k_ML2]=LeeCarterFit(D,E,a_ML,b_ML,k_ML);
  
%%%%%%%
[ MaxErr,x_MaxErr,t_MaxErr,Err2,L2,Rx,maxR2,x_maxR2,R2] = errorfkt( M,D,E,a_W,b_W,k_W );
