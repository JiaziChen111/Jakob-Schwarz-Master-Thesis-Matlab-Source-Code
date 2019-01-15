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

norming=0;


%load data with above parameters
BuildCentralDeathMatrix; 


  W=ones(x,t);
  index42=lnM==log(42);
  W(index42)=0;


BuildGDPvector;
BuildDrugData2;
alk=alk_2;
   
%alk_2  from 1950-2005 1:56
%rauch  from 1950-2003 1:54
%g_full from 1950-2016 1:67


t_future=50;
n_paths=10;
n=3;
g_future=zeros(n,t_future);
 
Md_g=arima(1,1,2);   
xxx=length(alk);
    EstMd_g = estimate(Md_g,alk(21:xxx)');
    h=forecast(EstMd_g,t_future+10,'Y0', alk(21:xxx)');
    H=h<0;
    h(H)=0;
    
     alky=figure;
     plot(linspace(1950,2005,56),alk,'g')
     hold on;
     plot(linspace(2005,2015,11),[alk(end),h(1:10)'],'r')
     hold on
  
    Md_g=arima(1,1,2);   

    EstMd_g = estimate(Md_g,alk');
    h2=forecast(EstMd_g,t_future+10,'Y0', alk');
    H=h2<0;
    h2(H)=0;
plot(linspace(2005,2015,11),[alk(end),h2(1:10)'],'b')
     hold on
     plot(linspace(2015,2015+t_future,t_future+1),h2(10:end),'b--')
     hold on
        plot(linspace(2015,2015+t_future,t_future+1),h(10:end),'r--')
  title('Daily alkohol consumption');
  legend('Historical','Forecast based on 1970-2005','Forecast based on 1950-2005')
 set(alky,'Units','Inches');
 pos = get(alky,'Position');
 set(alky,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 print(alky,'alky2','-dpdf','-r0')
    
alk=[alk_2(7:end),h(1:10)'];
g_future(2,:)=h(11:end);

if sex==3
    rauch=rauch_f;
    rauch_org=rauch_f_org;
elseif sex==4
    rauch=rauch_m;
    rauch_org=rauch_m_org;
elseif sex==5
    rauch=rauch_ges;
    rauch_org=rauch_ges_org;
end
    %rauch  from 1950-2003 1:54
    h=interp1(x_org,rauch_org',linspace(2004,2015+t_future,t_future+12),'linear','extrap');
    xxx=h(13:end);
    xxx=length(xxx);
    hdecrease=h(13:end)+linspace(1,xxx,xxx)*(-0.5);
    H=hdecrease<0;
    hdecrease(H)=0;
    
    hincrease=h(13:end)+linspace(1,xxx,xxx)*(0.5);
    H=hdecrease>100;
    hincrease(H)=100;
    
%     uu=figure;
%     plot(linspace(1950,2003,54),rauch,linspace(2004,2015+t_future,t_future+12),h,':.');
%     hold on
%     plot(linspace(2004,2015+t_future,t_future+12),h,'r:.');
%     plot(linspace(2016,2015+t_future,t_future),hdecrease,'b:.');
%     hold on
%     plot(linspace(2016,2015+t_future,t_future),hincrease,'b:.');
%    % axis([1950,2015+t_future 20 30])
%      title ('females smoking percentage')
%      set(uu,'Units','Inches');
%      pos = get(uu,'Position');
%      set(uu,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%      
%      print(uu,'uu2','-dpdf','-r0')
     rauch=[rauch(7:end)',h(1:12)];
     g_future(3,:)=h(13:end);
   
g=[g_t';rauch;alk];

     Md_g=arima(2,1,1);
     xxx=length(r_full);
     EstMd_g = estimate(Md_g,r_full(41:xxx));
     h=forecast(EstMd_g,t_future-1,'Y0', r_full(41:xxx));
     g_future(1,1)=(g_full(end));
     g_future(1,2)=(1+h(1))*exp(g_full(end));
     for i=3:t_future
         g_future(1,i)=(1+h(i-1))*g_future(1,i-1);
     end
         g_future(1,2:end)= log(g_future(1,2:end))';
    
%      h=interp1(linspace(1950,2016,67),g_full,linspace(2017,2015+t_future,t_future),'linear','extrap');
         
         
           gforecast=figure;
             plot(linspace(1950,1950+67,67),g_full(1:end),'g')
             hold on
             plot(linspace(2017,2017+t_future-1,t_future-1),g_future(1,2:end),'r')
%              
 Md_g=arima(2,1,1);
     
     EstMd_g = estimate(Md_g,r_full);
     h=forecast(EstMd_g,t_future-1,'Y0', r_full);
     g_future(1,1)=(g_full(end));
     g_future(1,2)=(1+h(1))*exp(g_full(end));
     for i=3:t_future
         g_future(1,i)=(1+h(i-1))*g_future(1,i-1);
     end
         g_future(1,2:end)= log(g_future(1,2:end))';



hold on
  plot(linspace(2017,2017+t_future-1,t_future-1),g_future(1,2:end),'b')
%         
 Md_g=arima(2,1,1);
     
     EstMd_g = estimate(Md_g,r_full(40:xxx));
     h=forecast(EstMd_g,t_future-1,'Y0', r_full(40:xxx));
     g_future(1,1)=(g_full(end));
     g_future(1,2)=(1+h(1))*exp(g_full(end));
     for i=3:t_future
         g_future(1,i)=(1+h(i-1))*g_future(1,i-1);
     end
         g_future(1,2:end)= log(g_future(1,2:end))';



hold on
  plot(linspace(2017,2017+t_future-1,t_future-1),g_future(1,2:end),'k')
%              plot(linspace(2017,2015+t_future,t_future),h,'b')
             title('log NI per capita')
             legend('Historical','Forecast based on 1950-2016','Forecast based on 1990-2016' ,'Forecast based on 1989-2016' ,'location','best');
            set(gforecast,'Units','Inches');
pos = get(gforecast,'Position');
set(gforecast,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gforecast,'gforecast4','-dpdf','-r0')

            