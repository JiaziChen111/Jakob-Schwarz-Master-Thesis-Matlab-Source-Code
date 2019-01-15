  %following code plots the error remaining in the reestimated k_LC
    
    tt=linspace(t_min,t_max,t);
    err_ML=zeros(1,t);
    for j=1:t
        for i=1:x
            err_ML(j)=err_ML(j)+E(i,j)*exp((a_ML(i)+b_ML(i)*k_ML(j)));
        end
        err_ML(j)=abs(sum (D(:,j))-err_ML(j)); %sum of all deaths year t-SumE);
    end
    err_wLC=zeros(1,t);
    for j=1:t
        for i=1:x
            err_wLC(j)=err_wLC(j)+E(i,j)*exp((a_wLC(i)+b_wLC(i)*k_wLC(j)));
        end
        err_wLC(j)=abs(sum (D(:,j))-err_wLC(j)); %sum of all deaths year t-SumE);
    end
     err_DwLC=zeros(1,t);
    for j=1:t
        for i=1:x
            err_DwLC(j)=err_DwLC(j)+E(i,j)*exp((a_DwLC(i)+b_DwLC(i)*k_DwLC(j)));
        end
        err_DwLC(j)=abs(sum (D(:,j))-err_DwLC(j)); %sum of all deaths year t-SumE);
    end
    xz=figure('name','Fits of Lee-Carter Modell second step');
  
  
    plot(tt,err_wLC,'r')
    hold on;
    plot(tt,err_DwLC,'k')
    hold on;
    plot(tt,err_ML,'b')
    title('err=|Sum_x D(x,t)-Sum_x N(x,t)*e^{(a_x+k_t*b_x)}|')
    legend('err_{1-0}','err_{death}','err_{ML}')
    ylabel('absolut error')
    
    set(xz,'Units','Inches');
pos = get(xz,'Position');
set(xz,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xz,'errorsML','-dpdf','-r0')

hhh=figure;
subplot(2,2,1);
plot(linspace(1,x,x),a_wLC2(1:x),'r');
hold on;
plot(linspace(1,x,x),a_DwLC(1:x),'k');
hold on;
plot(linspace(1,x,x),a_ML(1:x),'b');
title('a_x');

subplot(2,2,2);
plot(linspace(1,x,x),b_wLC2,'r');
hold on;
plot(linspace(1,x,x),b_DwLC,'k');
hold on;
plot(linspace(1,x,x),b_ML,'b');
title('b_x');

subplot(2,2,3:4);
plot(linspace(1,t,t),k_wLC2,'r');
hold on;
plot(linspace(1,t,t),k_DwLC,'k');
hold on;
plot(linspace(1,t,t),k_ML,'b');

title('k_t');
legend('1-0-weights','Death-weights','ML');
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh,'axbxkt_ML','-dpdf','-r0')

M_wLC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC(i,j)=exp(a_wLC(i)+b_wLC(i)*k_wLC(j));
        end
    end
    
M_D=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D(i,j)=exp(a_DwLC(i)+b_DwLC(i)*k_DwLC(j));
        end
    end    
  M_wLC2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC2(i,j)=exp(a_wLC2(i)+b_wLC2(i)*k_wLC2(j));
        end
    end
M_ML=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML(i,j)=exp(a_ML(i)+b_ML(i)*k_ML(j));
        end
    end
    
    lnM_wLC = log(M_wLC);
lnM_wLC2 = log(M_wLC2);
lnM_D = log(M_D);
lnM_ML = log(M_ML);

zz=figure('name','m(x,t) Fits of Lee-Carter Modell 1');
ax1 = subplot(4,1,1);
%plots development over time for a specific age
age=1; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_wLC2(age,:),'r')

hold on
plot(time,M_ML(age,:),'b')
h1=num2str(ages(age));


%plots development over time for a specific age
age=66; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_wLC2(age,:),'r')

hold on
plot(time,M_ML(age,:),'b')

h=num2str(ages(age));
title(['m(x,t) over time for age ' h1 ' (steeper) and ' h ])


ax12 = subplot(4,1,2);
%plots development over time for a specific age
age=41; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_wLC2(age,:),'r')

hold on
plot(time,M_ML(age,:),'b')
h1=num2str(ages(age));


%plots development over time for a specific age
age=21; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_wLC2(age,:),'r')

hold on
plot(time,M_ML(age,:),'b')
h=num2str(ages(age));
title(['m(x,t) over time for age ' h1 ' (upper) and ' h ])


ax15 = subplot(4,1,3);
%plots development over time for a specific age
age=101; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_wLC2(age,:),'r')

hold on
plot(time,M_ML(age,:),'b')
h=num2str(ages(age));
title(['m(x,t) over time for age ' h])

ax188 = subplot(4,1,4);
%plots development over time for a specific age
age=109; %1:x enstpricht x_min:x_max %%%%%%%%%%%%%%%%%%% choose
time=linspace(t_min,t_max,t);
ages=linspace(x_min,x_max,x);

index42=M(age,:)~=42;
tempM=M(age,:);

plot(time(index42),tempM(index42),'g')
hold on
plot(time,M_wLC2(age,:),'r')

hold on
plot(time,M_ML(age,:),'b')
h=num2str(ages(age));
title(['m(x,t) over time for age ' h])
legend('historical','1-0','ML')
set(zz,'Units','Inches');
pos = get(zz,'Position');
set(zz,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(zz,'mxt_ML','-dpdf','-r0')


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
plot(ages(view),lnM_wLC2(view,year),'r')
hold on
plot(ages(view),lnM_ML(view,year),'b')
hold on
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
plot(ages(view),lnM_wLC2(view,year4),'r')
hold on
plot(ages(view),lnM_ML(view,year4),'b')
hold on

h=num2str(year+t_min-1);
h1= num2str( ages(view_start));
h2= num2str(ages(view_end));
title(['ln[m(x,t)] for ages ' h1 ' to ' h2 ' in 1956 (upper) and 2015'])
legend('historical','1-0','ML')
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh,'mxt_ML2','-dpdf','-r0')
