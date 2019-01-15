clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=104; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=0; %if 1 then all values <= in M are replaced with 42

%load data with above parameters
BuildCentralDeathMatrix; 


%%%Lee Carter Model Fitting%%%
    
   
    [a_SVD,b_SVD,k_SVD]=SVDFit(lnM);    
    [a_LC,b_LC,k_LC]=LeeCarterFit(D,E,a_SVD,b_SVD,k_SVD);
    
    %following code plots the error remaining in the reestimated k_LC
    
    tt=linspace(t_min,t_max,t);
    err_LC=zeros(1,t);
    for j=1:t
        for i=1:x
            err_LC(j)=err_LC(j)+E(i,j)*exp((a_LC(i)+b_LC(i)*k_LC(j)));
        end
        err_LC(j)=abs(sum (D(:,j))-err_LC(j)); %sum of all deaths year t-SumE);
    end
    err_SVD=zeros(1,t);
    for j=1:t
        for i=1:x
            err_SVD(j)=err_SVD(j)+E(i,j)*exp((a_SVD(i)+b_SVD(i)*k_SVD(j)));
        end
        err_SVD(j)=abs(sum (D(:,j))-err_SVD(j)); %sum of all deaths year t-SumE);
    end
  
    
    max(err_LC) %0.0018 for males in west germany data 1958-2015 0-100+
    sum(err_LC) %.0142
    sqrt(sum(err_LC.^2)) % 0.0029
    
    n=size(M);
x=n(1);
t=n(2);
M_LC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_LC(i,j)=exp(a_LC(i)+b_LC(i)*k_LC(j));
        end
    end
    
  M_SVD=zeros(x,t);
    for i=1:x
        for j=1:t
            M_SVD(i,j)=exp(a_SVD(i)+b_SVD(i)*k_SVD(j));
        end
    end

lnM_SVD = log(M_SVD);
lnM_LC = log(M_LC);

figure('name','m(x,t) Fits of Lee-Carter Modell 1');
ax1 = subplot(3,1,1);
%plots development over time for a specific age
age=1; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_LC(age,:),'r')
hold on
plot(time,M_SVD(age,:),'k')
hold on
h1=num2str(ages(age));


%plots development over time for a specific age
age=66; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_LC(age,:),'r')
hold on
plot(time,M_SVD(age,:),'k')
hold on
h=num2str(ages(age));
title(['m(x,t) over time for age ' h1 ' (steeper) and ' h ])


ax12 = subplot(3,1,2);
%plots development over time for a specific age
age=41; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_LC(age,:),'r')
hold on
plot(time,M_SVD(age,:),'k')
hold on
h1=num2str(ages(age));


%plots development over time for a specific age
age=21; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_LC(age,:),'r')
hold on
plot(time,M_SVD(age,:),'k')
hold on
h=num2str(ages(age));
title(['m(x,t) over time for age ' h1 ' (upper) and ' h ])


ax15 = subplot(3,1,3);
%plots development over time for a specific age
age=x; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_LC(age,:),'r')
hold on
plot(time,M_SVD(age,:),'k')
hold on
h=num2str(ages(age));
title(['m(x,t) over time for age ' h])
legend('historical','LC','SVD')

hhh= figure('name','m(x,t) Fits of Lee-Carter Modell 2');

%plots development over age for a specific year
year=1; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
year1=15; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
year2=30; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
year3=45; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
year4=t; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
view_start=1; %1 for start age x_min %%%%%%%%%%%%%%%%%%%
view_end=x; %x for end age x-1, i.e. 100+ for x=101 %%%%%%%%%%%
view=linspace(view_start,view_end,view_end-view_start+1);

templnM=lnM(view,year);
indexln42=templnM~=log(42);

templnM1=lnM(view,year1);
indexln421=templnM1~=log(42);

templnM2=lnM(view,year2);
indexln422=templnM2~=log(42);

templnM3=lnM(view,year3);
indexln423=templnM3~=log(42);

templnM4=lnM(view,year4);
indexln424=templnM4~=log(42);

tempview=ages(view);

plot(tempview(indexln42),templnM(indexln42),'g')
hold on
plot(ages(view),lnM_LC(view,year),'r')
hold on
plot(ages(view),lnM_SVD(view,year),'k')
hold on
plot(linspace(1,x,x),a_LC,'b');
%plot(tempview(indexln421),templnM1(indexln421),'g')
%hold on
%plot(ages(view),lnM_LC(view,year1),'r')
%hold on
%plot(ages(view),lnM_SVD(view,year1),'k')
%hold on
%plot(tempview(indexln422),templnM2(indexln422),'g')
%hold on
%plot(ages(view),lnM_LC(view,year2),'r')
%hold on
%plot(ages(view),lnM_SVD(view,year2),'k')
%hold on
%plot(tempview(indexln423),templnM3(indexln423),'g')
%hold on
%plot(ages(view),lnM_LC(view,year3),'r')
%hold on
%plot(ages(view),lnM_SVD(view,year3),'k')
%hold on
plot(tempview(indexln424),templnM4(indexln424),'g')
hold on
plot(ages(view),lnM_LC(view,year4),'r')
hold on
plot(ages(view),lnM_SVD(view,year4),'k')
hold on

h=num2str(year+t_min-1);
h1= num2str( ages(view_start));
h2= num2str(ages(view_end));
title(['ln[m(x,t)] for ages ' h1 ' to ' h2 ' in 1956 (upper) and 2015'])
legend('historical','LC','SVD','a_x')
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh,'filename','-dpdf','-r0')



eps_LC=abs(M-M_LC);
eps_SVD=abs(M-M_SVD);
max(max(eps_LC) )
sum(sum(eps_LC) )/(x*t)
sqrt(sum(sum(eps_LC.^2)) )/(x*t)
max(max(eps_SVD) )
sum(sum(eps_SVD) )/(x*t)
sqrt(sum(sum(eps_SVD.^2)) )/(x*t)
%%%%%%%%%%
hhh=figure;
subplot(2,2,1);
plot(linspace(1,x,x),a_SVD,'k');
hold on;
plot(linspace(1,x,x),a_LC,'r');
hold on;
title('a_x');

subplot(2,2,2);
plot(linspace(1,x,x),b_SVD,'k');
hold on;
plot(linspace(1,x,x),b_LC,'r');
hold on;
title('b_x');

subplot(2,2,3:4);
plot(linspace(1,t,t),k_SVD,'k');
hold on;
plot(linspace(1,t,t),k_LC,'r');
hold on;

title('k_t');
legend('LC','SVD');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh,'axbxkt','-dpdf','-r0')

%clearvars h h1 h2 n age ages ax1 ax2 i j time view view_end view_start
%clearvars year index42 tempM templnM indexln42 tempview

errh=figure;
subplot(2,2,1)
Merr= M- M_LC ;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title(' ')

subplot(2,2,2)
Merr= M- M_LC ;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title(' ')

subplot(2,2,3)
Merr= M- M_LC ;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')
title(' ')

subplot(2,2,4)
Merr= M- M_LC ;
xscale = [1956 2015];
yscale = [0 104];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
title(' ')
suptitle('Fitting error colormap')
set(errh,'Units','Inches');
pos = get(errh,'Position');
set(errh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errh,'errh','-dpdf','-r0')

diff=abs(k_SVD-k_LC);

[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUx,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_SVD,b_SVD,k_SVD );
[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,FVUxLC,maxFVUx,x_maxFVUx,FVU ] = errorfkt( M,D,E,a_SVD,b_SVD,k_LC );
ErrMat=[ MaxErr,MinErr,x_MaxErr,t_MaxErr,x_MinErr,t_MinErr,MeanErr,Err2,L2,maxFVUx,x_maxFVUx,FVU ];

MeanFVUx=zeros(22,1);
MeanFVUx(1)=FVUx(1);
MeanFVUx(22)=sum(FVUx(102:105) )/4;
for i=1:20
    MeanFVUx(i+1)=sum(FVUx(2+5*(i-1):6+5*(i-1)) )/5;
end

xzf=figure('name','FVU_x');
    
   plot(linspace(0,x,x),FVUx,'k')
   hold on
   plot(linspace(0,x,x),FVUxLC,'r')
   
    title('FVU_x')
    legend('SVD','LC')
    ylabel('FVU_x')
    
    axis([0 110 0 100])
    xlabel('x')
    
    set(xzf,'Units','Inches');
pos = get(xzf,'Position');
set(xzf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xzf,'errors','-dpdf','-r0')
    
