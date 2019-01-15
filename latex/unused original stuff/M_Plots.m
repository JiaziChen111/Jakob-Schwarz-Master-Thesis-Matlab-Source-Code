n=size(M);
x=n(1);
t=n(2);
M_LC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_LC(i,j)=exp(a_LC(i)+b_LC(i)*k_LC(j));
        end
    end

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
    
M_ML=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML(i,j)=exp(a_ML(i)+b_ML(i)*k_ML(j));
        end
    end
   
M_SVD=zeros(x,t);
    for i=1:x
        for j=1:t
            M_SVD(i,j)=exp(a_SVD(i)+b_SVD(i)*k_SVD(j));
        end
    end

lnM_SVD = log(M_SVD);
lnM_ML = log(M_ML);
lnM_LC = log(M_LC);
lnM_wLC = log(M_wLC);
lnM_wLC2 = log(M_wLC2);

figure('name','m(x,t) Fits of Lee-Carter Modell');
ax1 = subplot(2,1,1);
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
plot(time,M_ML(age,:),'b')
hold on
plot(time,M_wLC(age,:),'c')
hold on
plot(time,M_wLC2(age,:),'m')
h=num2str(ages(age));
title(ax1,['m(x,t) over time for age ' h])
legend('historical','LC','SVD','ML','wLC','wLC2')

ax2 = subplot(2,1,2);
%plots development over age for a specific year
year=t; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%%%%% choose
view_start=1; %1 for start age x_min %%%%%%%%%%%%%%%%%%%
view_end=x; %x for end age x-1, i.e. 100+ for x=101 %%%%%%%%%%%
view=linspace(view_start,view_end,view_end-view_start+1);

templnM=lnM(view,year);
indexln42=templnM~=log(42);
tempview=ages(view);

plot(tempview(indexln42),templnM(indexln42),'g')
hold on
plot(ages(view),lnM_LC(view,year),'r')
hold on
plot(ages(view),lnM_SVD(view,year),'k')
hold on
plot(ages(view),lnM_ML(view,year),'b')
hold on
plot(ages(view),lnM_wLC(view,year),'c')
hold on
plot(ages(view),lnM_wLC2(view,year),'m')
h=num2str(year+t_min-1);
h1= num2str( ages(view_start));
h2= num2str(ages(view_end));
title(ax2,['ln[m(x,t)] for ages ' h1 ' to ' h2 ' in year ' h])

clearvars h h1 h2 n age ages ax1 ax2 i j time view view_end view_start
clearvars year index42 tempM templnM indexln42 tempview