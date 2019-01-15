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

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_W2,b_W2,k_W );
ErrMat_W2=[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,maxFVUx,x_maxFVUx,FVU ];
FVUxW=FVUx;
[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_wLC2,b_wLC2,k_wLC2 );
ErrMat_wLC2=[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,maxFVUx,x_maxFVUx,FVU ];
FVUxwLC2=FVUx;
[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_ML2,b_ML2,k_ML);
ErrMat_ML=[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,maxFVUx,x_maxFVUx,FVU ];
FVUxML=FVUx;
[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_D2,b_D2,k_D2 );
ErrMat_D=[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,maxFVUx,x_maxFVUx,FVU ];
FVUxD=FVUx;

xzf2=figure('name','FVU_x');
 plot(linspace(0,x,x),FVUxW,'c')
   hold on   
plot(linspace(0,x,x),FVUxwLC2,'r')
      hold on
      plot(linspace(0,x,x),FVUxD,'k')
         hold on
  
   plot(linspace(0,x,x),FVUxML,'b')
 
 
    title('FVU_x')
    legend('w_{x,t}=t+x_{max}-x','1-0','Death-weights','ML')
    ylabel('FVU_x')
    
   % axis([0 115 0 100])
    xlabel('x')
    
set(xzf2,'Units','Inches');
pos = get(xzf2,'Position');
set(xzf2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xzf2,'errors','-dpdf','-r0')