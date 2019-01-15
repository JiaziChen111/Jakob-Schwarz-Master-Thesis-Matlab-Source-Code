n_C=size(q);
x_C=n_C(1);
t_C=n_C(2);
q_C=zeros(x_C,t_C);
    for i=1:x_C
        for j=1:t_C
            h=A_1(j)+A_2(j)*((x_min_C+i-1)+(j-1));
            q_C(i,j)=exp(h)/(1+exp(h));
        end
    end

qC=figure('name','q(x,t) Fit with Cairns Model');
ax1 = subplot(4,1,1);
%plots development over time for a specific age
age=1; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);
ages=linspace(x_min_C,x_max_C,x_C);

plot(time,q(age,:),'g')
hold on
plot(time,q_C(age,:),'r')
hold on
h=num2str(ages(age));

title(ax1,['q(x,t) over time for fixed age ' h])


ax9 = subplot(4,1,2);
%plots development over time for a specific age
age=21; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);
ages=linspace(x_min_C,x_max_C,x_C);

plot(time,q(age,:),'g')
hold on
plot(time,q_C(age,:),'r')
hold on
h=num2str(ages(age));

title(ax9,['q(x,t) over time for fixed age ' h])


ax1 = subplot(4,1,3);
%plots development over time for a specific age
age=x_C; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);
ages=linspace(x_min_C,x_max_C,x_C);

plot(time,q(age,:),'g')
hold on
plot(time,q_C(age,:),'r')
hold on
h=num2str(ages(age));
title(ax1,['q(x,t) over time for fixed age ' h])


ax2 = subplot(4,1,4);
%plots development over age for a specific year
year=t_C; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%% choose
view_start=1; %1 for start age x_min %%%%%%%%%%%%%%%%%%
view_end=x_C; %x for end age x-1, i.e. 100+ for x=101 %%%%%%%%%%%%%
view=linspace(view_start,view_end,view_end-view_start+1);

plot(ages(view),(log(q(view,year))),'g')
hold on
plot(ages(view),(log(q_C(view,year))),'r')
year1=1;
plot(ages(view),(log(q(view,year1))),'g')
hold on
plot(ages(view),(log(q_C(view,year1))),'r')

h=num2str(year+t_min_C-1);
h1x=num2str(year1+t_min_C-1);

h1= num2str( ages(view_start));
h2= num2str(ages(view_end));
title(ax2,['ln q(x,t) for ages ' h1 ' to ' h2 ' in year ' h1x ' and ' h ' (lower)'])
set(qC,'Units','Inches');
pos = get(qC,'Position');
set(qC,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(qC,'qC2','-dpdf','-r0')