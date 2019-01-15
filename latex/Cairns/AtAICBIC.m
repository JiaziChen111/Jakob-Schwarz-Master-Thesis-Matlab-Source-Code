clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =4; % chose sex 3 4 5 for female male total
x_max=110; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42

%load data with above parameters
BuildCentralDeathMatrix;

%%%Cairns model
   x_max_C=109; %male goes only to 109
   x_min_C= 70; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1, A_2,mu,C]=q_RegressionFit(q,x_min_C,x_max_C);

    AICBIC=zeros(6,2);
    Dt=[diff(A_1),diff(A_2)];
    for i=0:5
        Md=varm(2,i);
        EstMd = estimate(Md,Dt);
        results = summarize(EstMd);
        AICBIC(i+1,1) = results.AIC;
        AICBIC(i+1,2) = results.BIC;
    end
   fprintf([repmat('%.2f\t', 1, size(AICBIC, 2)) '\n'], AICBIC')

    