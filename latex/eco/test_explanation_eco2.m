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

BuildGDPvector;
BuildDrugData;

  W=ones(x,t);
  index42=lnM==log(42);
  W(index42)=0;

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
norming=0;

a_start=zeros(x,1);
b_start=ones(x,1);
k_start=ones(t,1);
y_start=zeros(x,n);

a_startML=zeros(x,1);
b_startML=ones(x,1);
k_startML=zeros(t,1);
y_startML=zeros(x,n);

[a_wLC,b_wLC,k_wLC,fval_wLC,exitflag_wLC,steps_wLC]=Fit_wLC(lnM,W);
[a_wLC_re,b_wLC_re,k_wLC_re,y_wLC_re]=explanation_renorm(a_wLC,b_wLC,k_wLC,g ,norming);
tic
[a_2step,b_2step,k_2step,y_2step,fval_2step,exitflag_2step,steps_2steps]=FitLC_Explanation_2step(lnM,g,norming,W);
time2step=toc;
[a_2stepD,b_2stepD,k_2stepD,y_2stepD,fval_2stepD,exitflag_2stepD,steps_2stepsD]=FitLC_Explanation_2step(lnM,g,norming,D);

tic
[a_ex_newton,b_ex_newton,k_ex_newton,y_ex_newton,fval_ex_newton,exitflag_ex_newton]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,W,norming);
timeLC_newton=toc;
[a_D,b_D,k_D,y_D,fval_D,exitflag_D]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,D,norming);
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

  [a_ML_newton2step,b_ML_newton2step,k_ML_newton2step,y_ML_newton2step] = LeeCarterFit_Explanation(D, E, a_ML_newton, b_ML_newton, k_ML_newton, y_ML_newton, g, norming);
 y0=reshape(y_ML_newton2step,x*n,1);
        x0 = [a_ML_newton2step;b_ML_newton2step;k_ML_newton2step;y0];
        fval_ML_newton2step=FUN(x0);
        
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

[a_2step2,b_2step2,k_2step2,y_2step2] = LeeCarterFit_Explanation(D, E, a_2step, b_2step, k_2step, y_2step, g, norming);
 y0=reshape(y_2step2,x*n,1);
        x0 = [a_2step2;b_2step2;k_2step2;y0];
        fval_2step2=FUN(x0);

norming=1;
[a_2stepN,b_2stepN,k_2stepN,y_2stepN,fval_2stepN,exitflag_2stepN,steps_2stepsN]=FitLC_Explanation_2step(lnM,g,norming,W);


lnM_2step=zeros(x,t);   
    for i=1:x
        for j=1:t
            lnM_2step(i,j)=(a_2step(i)+b_2step(i)*k_2step(j)+y_2step(i,:)*g(:,j));

        end
    end
  
     lnM_wLC=zeros(x,t);   
    for i=1:x
        for j=1:t
            lnM_wLC(i,j)=(a_wLC(i)+b_wLC(i)*k_wLC(j));

        end
    end
    lnM_wLC_newton=zeros(x,t);   
    for i=1:x
        for j=1:t
           lnM_wLC_newton(i,j)=(a_ex_newton(i)+b_ex_newton(i)*k_ex_newton(j)+y_ex_newton(i,:)*g(:,j));

        end
    end
     lnM_wLC_newtonD=zeros(x,t);   
    for i=1:x
        for j=1:t
           lnM_wLC_newtonD(i,j)=(a_D(i)+b_D(i)*k_D(j)+y_D(i,:)*g(:,j));

        end
    end
    lnM_wLC_newton3step=zeros(x,t);   
    for i=1:x
        for j=1:t
           lnM_wLC_newton3step(i,j)=(a_ex_newton2(i)+b_ex_newton2(i)*k_ex_newton2(j)+y_ex_newton2(i,:)*g(:,j));

        end
    end 
     lnM_ML_normalen=zeros(x,t);   
    for i=1:x
        for j=1:t
           lnM_ML_normalen(i,j)=(a_ML(i)+b_ML(i)*k_ML(j)+y_ML(i,:)*g(:,j));

        end
    end
     lnM_ML_newton=zeros(x,t);   
    for i=1:x
        for j=1:t
          lnM_ML_newton(i,j)=(a_ML_newton(i)+b_ML_newton(i)*k_ML_newton(j)+y_ML_newton(i,:)*g(:,j));

        end
    end
        lnM_ML_newton3step=zeros(x,t);   
    for i=1:x
        for j=1:t
          lnM_ML_newton3step(i,j)=(a_ML_newton2(i)+b_ML_newton2(i)*k_ML_newton2(j)+y_ML_newton2(i,:)*g(:,j));

        end
    end
     lnM_ML_normalen3step=zeros(x,t);   
    for i=1:x
        for j=1:t
           lnM_ML_normalen3step(i,j)=(a_ML2(i)+b_ML2(i)*k_ML2(j)+y_ML2(i,:)*g(:,j));

        end
    end
 M_2step=exp(lnM_2step);
 M_wLC=exp(lnM_wLC) ; 
 M_wLC_newton=exp(lnM_wLC_newton) ; 
 M_wLC_newtonD=exp(lnM_wLC_newtonD) ; 
 M_wLC_newton3step=exp(lnM_wLC_newton3step) ; 
 M_ML_normalen=exp(lnM_ML_normalen) ; 
 M_ML_newton=exp(lnM_ML_newton) ; 
 M_ML_newton3step=exp(lnM_ML_newton3step) ; 
 M_ML_normalen3step=exp(lnM_ML_normalen3step) ; 
%%%%%%%%%%
x_cut=x   ;   
hhhx=figure;

subplot(2,2,1);
plot(linspace(1-1,x_cut-1 ,x_cut ),a_2step(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut-1 ,x_cut ),a_wLC_re(1:x_cut ),'r');
hold on;
plot(linspace(1,x_cut-1 ,x_cut ),a_D(1:x_cut ),'k');

title('a_x');

subplot(2,2,2);
plot(linspace(1-1,x_cut-1 ,x_cut ),b_2step(1:x_cut ),'b');
hold on;
plot(linspace(1-1,x_cut-1 ,x_cut ),b_wLC_re(1:x_cut ),'r');
hold on;
plot(linspace(1-1,x_cut-1 ,x_cut ),b_D(1:x_cut ),'k');
title('b_x');

subplot(2,2,3);
plot(linspace(1,t ,t ),k_2step,'b');
hold on;
plot(linspace(1,t ,t ),k_wLC_re,'r');
hold on;
plot(linspace(1,t ,t ),k_D,'k');
title('k_t');

subplot(2,2,4);
plot(linspace(1-1,x_cut-1 ,x_cut ),y_2step(1:x_cut ),'b');
hold on;
plot(linspace(1-1,x_cut-1 ,x_cut ),y_wLC_re(1:x_cut ),'r');
hold on;
plot(linspace(1-1,x_cut-1 ,x_cut ),y_D(1:x_cut ),'k');
title('y_x');

legend({'3step','wLC-renormed','death-weights'},'location','best');
set(hhhx,'Units','Inches');
pos = get(hhhx,'Position');
set(hhhx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhx,'renorm','-dpdf','-r0')

%%%%%%%%%
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
print(hhh,'gdpabk2','-dpdf','-r0')
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
print(hhhhh,'gdpy2','-dpdf','-r0')

x_cut=x   ;   
xxx=figure;
  
aX_2step=zeros(x,1);
for i=1:x
   aX_2step(i)=sum(lnM_2step(i,:))/t;    
end
aX_2stepN=zeros(x,1);
for i=1:x
   aX_2stepN(i)=sum(lnM_2stepN(i,:))/t;    
end

subplot(2,2,1);
plot(linspace(1,x_cut ,x_cut ),a_wLC(1:x_cut ),'r');
hold on;
plot(linspace(1,x_cut ,x_cut ),aX_2step(1:x_cut ),'g');
hold on;
plot(linspace(1,x_cut ,x_cut ),aX_2stepN(1:x_cut ),'b');
hold on;


title('a^{*}_x');

subplot(2,2,2);

plot(linspace(1,x_cut ,x_cut ),b_wLC(1:x_cut ),'r');
hold on;
plot(linspace(1,x_cut ,x_cut ),b_2step(1:x_cut ),'g');
hold on;
plot(linspace(1,x_cut ,x_cut ),b_2stepN(1:x_cut ),'b');
hold on;

%plot(linspace(1,x_cut ,x_cut ),b_normalen(1:x_cut ),'g');


title('b_x');


subplot(2,2,3:4);
plot(linspace(1,t,t),k_wLC(1:t),'r');
hold on;
plot(linspace(1,t,t),k_2step(1:t),'g');
hold on;
plot(linspace(1,t,t),k_2stepN(1:t),'b');


%plot(linspace(1,t,t),k_normalen(1:t),'g');

title('k_t');
legend({'classic-wLC','3step','3step_{alt}'},'location','best');
set(xxx,'Units','Inches');
pos = get(xxx,'Position');
set(xxx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xxx,'gdpabkwLC','-dpdf','-r0')

errq=figure;

Merr= abs(M-M_2step);
MerrwLC=abs( M-M_wLC );
Merr(index42)=0;
MerrwLC(index42)=0;

subplot(1,2,1)

xscale = [1956 2015];
yscale = [0 110];
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([0 0.1]);
xlabel('year t')
ylabel('age x')
title('wLC')

subplot(1,2,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([0 0.1]);
xlabel('year t')
ylabel('age x')
title('2step')

suptitle('Absolut fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errGDPwLC','-dpdf','-r0')

errq=figure;
Merr= M_2step-M;
MerrwLC= M_wLC-M ;
Merr(index42)=0;
MerrwLC(index42)=0;

subplot(2,2,1)

xscale = [1956 2015];
yscale = [0 110];
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([-0.1 0.1]);
xlabel('year t')
ylabel('age x')
title('wLC');
subplot(2,2,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.1 0.1]);
xlabel('year t')
ylabel('age x')
title('2step');
subplot(2,2,3)
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')

suptitle('fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errGDP2step','-dpdf','-r0')

%%%%% relative

Merr= (M-M_2step)./M;
MerrwLC= (M-M_wLC)./M ;
Merr(index42)=0;
MerrwLC(index42)=0;

errq=figure;
subplot(2,2,1)

xscale = [1956 2015];
yscale = [0 110];
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([-0.2 0.2]);
xlabel('year t')
ylabel('age x')
title('wLC');
subplot(2,2,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.2 0.2]);
xlabel('year t')
ylabel('age x')
title('2step');
subplot(2,2,3)
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([-0.05 0.05]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.05]);
xlabel('year t')
ylabel('age x')

suptitle('relative fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'relerrGDP2step','-dpdf','-r0')

%%%%% percenatge

Merr= (M_2step./M)-1;
MerrwLC= (M_wLC./M)-1 ;
Merr(index42)=0;
MerrwLC(index42)=0;

errq=figure;

subplot(2,2,1)
xscale = [1956 2015];
yscale = [0 110];
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([-0.4 0.4]);
xlabel('year t')
ylabel('age x')
title('wLC');

subplot(2,2,2)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.4 0.4]);
xlabel('year t')
ylabel('age x')
title('2step');

subplot(2,2,3)
xscale = [1956 2015];
yscale = [0 110];
imagesc(xscale,yscale,MerrwLC)
colorbar
caxis([-0.2 0.2]);
xlabel('year t')
ylabel('age x')
title('wLC');

subplot(2,2,4)
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.2 0.2]);
xlabel('year t')
ylabel('age x')
title('2step');

suptitle('relative fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'percerrGDP2step','-dpdf','-r0')
%%%%%%%%%%%%%%
err_2step=zeros(1,t);
    for j=1:t
        for i=1:x
            for k=1:n
            err_2step(j)=err_2step(j)+E(i,j)*exp(a_2step(i)+b_2step(i)*k_2step(j)+y_2step(i,k)*g(k,j));
            end
        end       
        err_2step(j)=abs(sum (D(:,j))-err_2step(j)); %sum of all deaths year t-SumE);
    end
 err_ML=zeros(1,t);
    for j=1:t
        for i=1:x
            for k=1:n
            err_ML(j)=err_ML(j)+E(i,j)*exp(a_ML(i)+b_ML(i)*k_ML(j)+y_ML(i,k)*g(k,j));
            end
        end       
        err_ML(j)=abs(sum (D(:,j))-err_ML(j)); %sum of all deaths year t-SumE);
    end
  err_2stepD=zeros(1,t);
    for j=1:t
        for i=1:x
            for k=1:n
            err_2stepD(j)=err_2stepD(j)+E(i,j)*exp(a_2stepD(i)+b_2stepD(i)*k_2stepD(j)+y_2stepD(i,k)*g(k,j));
            end
        end       
        err_2stepD(j)=abs(sum (D(:,j))-err_2stepD(j)); %sum of all deaths year t-SumE);
    end
    err_D=zeros(1,t);
    for j=1:t
        for i=1:x
            for k=1:n
            err_D(j)=err_D(j)+E(i,j)*exp(a_D(i)+b_D(i)*k_D(j)+y_D(i,k)*g(k,j));
            end
        end       
        err_D(j)=abs(sum (D(:,j))-err_D(j)); %sum of all deaths year t-SumE);
    end
   
  xz=figure;
 
  tt=linspace(1,t,t);
    plot(tt,err_2step,'r')
    hold on;
   plot(tt,err_D,'b')
   hold on;
   plot(tt,err_ML,'k')
    title('err_t=|Sum_x D(x,t)-Sum_x N(x,t)*e^{(a_x+k_t*b_x+y_x*g_t)}|')
    legend('err_{2step}','err_{D}','err_{ML}','location','best')
    ylabel('absolut error')
    
    
    set(xz,'Units','Inches');
pos = get(xz,'Position');
set(xz,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xz,'GDPerrors','-dpdf','-r0')

 
ErrM=zeros(8,2);
[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_1x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_2step,D,E );
ErrM(2,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_2x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_wLC_newton,D,E );
ErrM(1,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_3x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_wLC_newtonD,D,E );
ErrM(3,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_4x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_wLC_newton3step,D,E);
ErrM(4,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_5x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_ML_normalen,D,E  );
ErrM(5,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_6x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_ML_newton,D,E );
ErrM(6,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_7x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_ML_newton3step,D,E);
ErrM(7,:)=[ Err2,R2];

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_8x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_ML_normalen3step,D,E );
ErrM(8,:)=[ Err2,R2];

fprintf([repmat('%.5f\t', 1, size(ErrM, 2)) '\n'], ErrM')

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,R_0x,maxFVUx,x_maxFVUx,R2 ] = errorfkt2( M,M_wLC,D,E );

xzf=figure;
subplot(2,1,1)
 plot(linspace(0,90,91),R_0x(1:91),'g')
  hold on; 
plot(linspace(0,90,91),R_1x(1:91),'r')
  hold on;
 % plot(linspace(0,90,91),R_2x(1:91),'k')
  hold on;
  %plot(linspace(0,110,x),R_3x,'k')
  hold on;
  plot(linspace(0,90,91),R_4x(1:91),'k')
  hold on;
  plot(linspace(0,90,91),R_5x(1:91),'b')
  hold on;
 % plot(linspace(0,110,x),R_6x,'c')
  hold on;
  %plot(linspace(0,110,x),R_7x,'m')
  hold on;
 % plot(linspace(0,110,x),R_8x,'w')
   %hold on
   %plot(linspace(100,110,11),R_wLC,'r')
   %hold on
   %plot(linspace(100,110,11),R_D,'k')
   %hold on
   %plot(linspace(100,110,11),R_ML,'b')
   
   %legend('3step','wLC\_newton','	wLC\_newton\_D	','	wLC\_newton\_3step','	ML\_normalen',	'ML\_newton'	,	'ML\_normalen\_3step'	,	'ML\_newton\_3step','location','best')
   legend('classic 1-0-wLC','3step','wLC\_newton\_3step','ML\_normalen'	,'location','best')
   ylabel('FVU_x')
    xlabel('x')
subplot(2,1,2)    
 plot(linspace(90,110,110-90+1),R_0x(91:end),'g')
  hold on;  
plot(linspace(90,110,110-90+1),R_1x(91:end),'r')
  hold on;
%  plot(linspace(90,110,110-90+1),R_2x(91:end),'k')
  hold on;
  %plot(linspace(0,110,x),R_3x,'k')
  hold on;
  plot(linspace(90,110,110-90+1),R_4x(91:end),'k')
  hold on;
  plot(linspace(90,110,110-90+1),R_5x(91:end),'b')
  hold on;
 % plot(linspace(0,110,x),R_6x,'c')
  hold on;
  %plot(linspace(0,110,x),R_7x,'m')
  hold on;
 % plot(linspace(0,110,x),R_8x,'w')
   %hold on
   %plot(linspace(100,110,11),R_wLC,'r')
   %hold on
   %plot(linspace(100,110,11),R_D,'k')
   %hold on
   %plot(linspace(100,110,11),R_ML,'b')
   
   %legend('3step','wLC\_newton','	wLC\_newton\_D	','	wLC\_newton\_3step','	ML\_normalen',	'ML\_newton'	,	'ML\_normalen\_3step'	,	'ML\_newton\_3step','location','best')
    ylabel('FVU_x')
    xlabel('x')
   suptitle('Fractional Variance Unexplained')
%title('FVU_x')
 
    set(xzf,'Units','Inches');
pos = get(xzf,'Position');
set(xzf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xzf,'FVU_gdp','-dpdf','-r0')
 
AICBIC=zeros(6,2);
h=diff(k_2stepN)';
    Dt=[h(2:end);diff(diff(g))]';
    for i=0:5
        Md=varm(2,i);
        EstMd = estimate(Md,Dt);
        results = summarize(EstMd);
        AICBIC(i+1,1) = results.AIC;
        AICBIC(i+1,2) = results.BIC;
    end
   fprintf([repmat('%.2f\t', 1, size(AICBIC, 2)) '\n'], AICBIC')

 AICBIC_kt=zeros(9,2);
 k=1 ;  
    for i=0:2
        for j=0:2
        
            Md=arima(i,1,j);
        EstMd = estimate(Md,r_full);
        results = summarize(EstMd);
        AICBIC_kt(k,1) = results.AIC;
        AICBIC_kt(k,2) = results.BIC;
        k=k+1;
        end
    end
   fprintf([repmat('%.2f\t', 1, size(AICBIC_kt, 2)) '\n'], AICBIC_kt')

   
dk=figure('name','Identifying d');
subplot(2,2,1);
h=length(k_wLC)-1;
plot(linspace(1,h,h),diff(g),'r');
title('dg_t');

subplot(2,2,2) 
autocorr(diff(g))

ylabel(' ')
title('Sample Autocorr.Function for dg_t')

subplot(2,2,3) 
h=length(diff(diff(g)));
plot(linspace(1,h,h),diff(diff(g)),'r');

title('dd g_{t}');
subplot(2,2,4)
autocorr(diff(diff(g)))

ylabel(' ')
title('d^2 g_{t}')

subplot(2,2,4) 
title('Sample Autocorr.Function for ddg_t')

set(dk,'Units','Inches');
pos = get(dk,'Position');
set(dk,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(dk,'dkgdp','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%
dk=figure('name','Identifying d');
subplot(2,2,1);
autocorr(diff(k_2step))
title('dk\_3step')

subplot(2,2,2) 

parcorr(diff(k_2step))

title('dk\_3step_{alt}')
subplot(2,2,3) 

autocorr(diff(k_2stepN))
title(' ')

subplot(2,2,4)

parcorr(diff(k_2stepN))

%%%%%%%%%%%%%%%%%%
dk=figure('name','Identifying d');
subplot(2,2,1);
plot(linspace(1951,2016,66),r_full);
title('r_t')

subplot(2,2,3) 

autocorr(r_full)
title(' ')

subplot(2,2,2);
plot(linspace(1952,2016,65),diff(r_full));
title('dr_t')

subplot(2,2,4) 

autocorr(diff(r_full))
title(' ')
set(dk,'Units','Inches');
pos = get(dk,'Position');
set(dk,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(dk,'r_t','-dpdf','-r0')


dk=figure('name','Identifying d');
subplot(2,1,1);
autocorr(diff(diff(g)))
title('dd g_t')

subplot(2,1,2) 

parcorr(diff(diff(g)))

title(' ')


set(dk,'Units','Inches');
pos = get(dk,'Position');
set(dk,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(dk,'dkgdp4','-dpdf','-r0')