n_C=size(q);
x_C=n_C(1);
t_C=n_C(2);
q_C=zeros(x_C,t_C);
q_C_future=zeros(x_C,t_future_C);
    for i=1:x_C
        for j=1:t_C+t_future_C
            if j<=t_C
                h=A_1(j)+A_2(j)*((x_min_C+i-1)+(j-1));
                q_C(i,j)=exp(h)/(1+exp(h));
            else            
                h=A_future_mean(j-t_C,1)+A_future_mean(j-t_C,2)*((x_min_C+i-1)+(j-1));
                q_C_future(i,j-t_C)=exp(h)/(1+exp(h));
            end
        end
    end
lnq_C=log(q_C);
lnq_full=log(q_full);
lnq_C_future=log(q_C_future);

%calc alpha/2 and 1-alpha/2 quantil for q_C from sample path
    alpha=0.1;
    q_C_future_UB=zeros(x_C,t_future_C);
    q_C_future_LB=zeros(x_C,t_future_C);
    H=zeros(n_paths_C,1);
    for i=1:x_C
        for j=1:t_future_C
            for r=1:n_paths_C
                h=A_future{r}(j,1)+A_future{r}(j,2)*((x_min_C+i-1)+(j-1+t_C));
                H(r,1)=exp(h)/(1+exp(h));
            end
       q_C_future_LB(i,j) = prctile(H',100*(alpha/2),2);
       q_C_future_UB(i,j) = prctile(H',100*(1-alpha/2),2);    
        end
    end

    lnq_C_future_LB=log(q_C_future_LB);
    lnq_C_future_UB=log(q_C_future_UB);


fig1=figure('name','q(x,t) Fit and Forecast with Cairns Model');
%plots development over time for a specific age
age=x_max_C; %%%%%%%%%%%%%%choose
time_full=linspace(t_start_C,t_end_C,t_full_C);
time=linspace(t_min_C,t_max_C,t_C);
time_future=linspace(t_max_C,t_max_C+t_future_C,t_future_C+1);
ages=x_min_C:x_max_C;

plot(time_full,q_full(age+1,:),'g')
hold on
plot(time,q_C(age+1-x_min_C,:),'r')
hold on
h=num2str(age);
title(['q(x,t) over time for fixed age ' h])
legend('historical','Cairns')
plot(time_future,[q_C(end),q_C_future_UB(age+1-x_min_C,:,:)],'r:')
plot(time_future,[q_C(end),q_C_future_LB(age+1-x_min_C,:,:)],'r:')
plot(time_future,[q_C(end),q_C_future(age+1-x_min_C,:,:)],'r--')

fig2=figure;
%plots development over age for a specific forecasted year
year=t_future_C; %beetween 1:t_future_C %%%%%%%%%%% choose
if year+t_C<=t_full_C
    plot(ages,q_full(ages,year+t_C),'g')
    hold on
end
plot(ages,q_C_future(ages-x_min_C+1,year),'r--')
hold on
plot(ages,q_C_future_UB(ages-x_min_C+1,year),'r:')
hold on
plot(ages,q_C_future_LB(ages-x_min_C+1,year),'r:')
h=num2str(year+t_max_C);
h1= num2str(ages(1));
h2= num2str(ages(end));
title(['ln q(x,t) for ages ' h1 ' to ' h2 ' in fixed year ' h])

clearvars h h1 h2 n_C age ages ax1 ax2 i j time view view_end view_start year H h r 
clearvars q_C_future_UB q_C_future_LB lnq_C lnq_full lnq_C_future_LB lnq_C_future_UB