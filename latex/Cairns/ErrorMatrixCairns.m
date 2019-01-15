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

%%%Cairns model
   x_max_C=110; %for a new maximum age group, i .e. 100 instead of 110
   x_min_C= 40; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1, A_2,mu,C]=q_RegressionFit(q,x_min_C,x_max_C);

q_C=zeros(x_C,t_C);
    for i=1:x_C
        for j=1:t_C
            h=A_1(j)+A_2(j)*((x_min_C+i-1)+(j-1));
            q_C(i,j)=exp(h)/(1+exp(h));
        end
    end




%no weights but second step adjusted
%use implied M by q
lnM=log(-log(1-q_full));
[a_wLC,b_wLC,k_wLC,steps_wLC,fval_wLC]=Fit_wLC(lnM);
[a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);

%weights given by D
[a_D,b_D,k_D,steps_D,fval_D]=Fit_wLC(lnM,D);
[a_D2,b_D2,k_D2]=LeeCarterFit(D,E,a_D,b_D,k_D);

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
M_wLC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC(i,j)=exp(a_wLC(i)+b_wLC(i)*k_wLC(j));
        end
    end
M_wLC2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC2(i,j)=exp(a_wLC2(i)+b_wLC2(i)*k_wLC2(j));
        end
    end    
M_D=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D(i,j)=exp(a_D(i)+b_D(i)*k_D(j));
        end
    end    
    M_D2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D2(i,j)=exp(a_D2(i)+b_D2(i)*k_D2(j));
        end
    end    
  M_W=zeros(x,t);
    for i=1:x
        for j=1:t
            M_W(i,j)=exp(a_W(i)+b_W(i)*k_W(j));
        end
    end
M_W2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_W2(i,j)=exp(a_W2(i)+b_W2(i)*k_W2(j));
        end
    end
M_ML=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML(i,j)=exp(a_ML(i)+b_ML(i)*k_ML(j));
        end
    end
M_ML2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML2(i,j)=exp(a_ML2(i)+b_ML2(i)*k_ML2(j));
        end
    end
    
q_wLC=1-exp(-M_wLC);
q_wLC2=1-exp(-M_wLC2);
q_W=1-exp(-M_W);
q_W2=1-exp(-M_W2);
q_D=1-exp(-M_D);
q_D2=1-exp(-M_D2);
q_ML=1-exp(-M_ML);
q_ML2=1-exp(-M_ML2);
%q_M=1-exp(-M);

%[ MaxErr,x_MaxErr,t_MaxErr,Err2,L2,Rx,maxR2,x_maxR2,R2] = errorfkt( M,D,E,a_W,b_W,k_W );


ErrM=zeros(9,5);
[ MaxErr,MinErr,MeanErr,Err2,R2,R_C]=errorfkt_q( q,q_C );
ErrM(1,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_wLC]=errorfkt_q( q,q_wLC(x_min_C+1:end,:) );
ErrM(2,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_wLC2]=errorfkt_q( q,q_wLC2(x_min_C+1:end,:) );
ErrM(3,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_D]=errorfkt_q( q,q_D(x_min_C+1:end,:) );
ErrM(4,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_D2]=errorfkt_q( q,q_D2(x_min_C+1:end,:) );
ErrM(5,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_W]=errorfkt_q( q,q_W(x_min_C+1:end,:) );
ErrM(6,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_W2]=errorfkt_q( q,q_W2(x_min_C+1:end,:) );
ErrM(7,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_ML]=errorfkt_q(q,q_ML(x_min_C+1:end,:) );
ErrM(8,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_ML2]=errorfkt_q( q,q_ML2(x_min_C+1:end,:) );
ErrM(9,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

fprintf([repmat('%.5f\t', 1, size(ErrM, 2)) '\n'], ErrM')

xzf=figure('name','FVU_x');
    subplot(2,1,1)
   plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_C,'g')
   hold on
   plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_wLC,'r')
   hold on
    plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_D,'k')
   hold on
%  plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_W,'c')
 
    ylabel('FVU_x')
    xlabel('x')
   
  legend({'Cairns_{40-110}','1-0-weights','Death-weights,','W'},'Location','northwest')  
   
 %%%%%%%%%%%%%
  x_max_C=110; %for a new maximum age group, i .e. 100 instead of 110
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

q_C=zeros(x_C,t_C);
    for i=1:x_C
        for j=1:t_C
            h=A_1(j)+A_2(j)*((x_min_C+i-1)+(j-1));
            q_C(i,j)=exp(h)/(1+exp(h));
        end
    end


%[ MaxErr,x_MaxErr,t_MaxErr,Err2,L2,Rx,maxR2,x_maxR2,R2] = errorfkt( M,D,E,a_W,b_W,k_W );

[ MaxErr,MinErr,MeanErr,Err2,R2,R_C]=errorfkt_q( q,q_C );
 %%%%%%%%%%%%
 
subplot(2,1,2)
  plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_C,'g')
   %hold on
   %plot(linspace(100,110,11),R_wLC,'r')
   %hold on
   %plot(linspace(100,110,11),R_D,'k')
   %hold on
   %plot(linspace(100,110,11),R_ML,'b')
   
   legend('Cairns_{70-110}')
   ylabel('FVU_x')
    xlabel('x')
    
   suptitle('FVU_x')
%title('FVU_x')
 
    set(xzf,'Units','Inches');
pos = get(xzf,'Position');
set(xzf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xzf,'FVU_q2','-dpdf','-r0')
    

