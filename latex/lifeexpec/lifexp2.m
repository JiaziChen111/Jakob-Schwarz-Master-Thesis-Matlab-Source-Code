clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=100; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
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


t_future=110;
n_paths=10;
n=3;
g_future=zeros(n,t_future);
g_future_paths=cell(n_paths,1);

Md_g=arima(1,1,2);   
    %mean path
    xxx=length(alk);
    EstMd_g = estimate(Md_g,alk(21:xxx)');
    h=forecast(EstMd_g,t_future+10,'Y0', alk(21:xxx)');
    H=h<0;
    h(H)=0;
    
    %paths
    res = infer(EstMd_g,alk');
    Ysim = simulate(EstMd_g,t_future+10,'NumPaths',n_paths,'Y0',alk','E0',res);
    H=Ysim<0;
    Ysim(H)=0; 
    for i=1:n_paths
      g_future_paths{i}=zeros(n,t_future);
      g_future_paths{i}(2,:)=Ysim(11:end,i);
    end
        
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
    H=h<0;
    h(H)=0;
    H=h>100;
    h(H)=100; 
    
    for i=1:n_paths
    xxx=length(h(13:end));    
    xrate=rand-0.5   ;
  
    hrand=h(13:end)+linspace(1,xxx,xxx)*(xrate);
    H=hrand<0;
    hrand(H)=0;
    H=hrand>100;
    hrand(H)=100;
    
    g_future_paths{i}(3,:)=hrand;
    
    end
  
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
     EstMd_g = estimate(Md_g,r_full(41:66));
     h=forecast(EstMd_g,t_future-1,'Y0', r_full(41:end));
     g_future(1,1)=g_full(end);
     g_future(1,2)=(1+h(1))*exp(g_full(end));
     for i=3:t_future
         g_future(1,i)=(1+h(i-1))*g_future(1,i-1);
     end
         g_future(1,2:end)= log(g_future(1,2:end))';
          
    %paths
    res = infer(EstMd_g, r_full(41:end));
    Ysim = simulate(EstMd_g,t_future-1,'NumPaths',n_paths,'Y0', r_full(41:end),'E0',res);
    for i=1:n_paths
      g_future_paths{i}(1,1)=g_full(end);
      g_future_paths{i}(1,2)=(1+Ysim(1,i))*exp(g_full(end));
      for ii=3:t_future
         g_future_paths{i}(1,ii)=(1+Ysim(ii-1,i))*g_future_paths{i}(1,ii-1);
      end
         g_future_paths{i}(1,2:end)= log(g_future_paths{i}(1,2:end))';
    end
 
a_start=zeros(x,1);
b_start=ones(x,1);
k_start=ones(t,1);
y_start=zeros(x,n);

a_startML=zeros(x,1);
b_startML=ones(x,1);
k_startML=zeros(t,1);
y_startML=zeros(x,n);

    

 a_x=cell(8,1);
 b_x=cell(8,1);
 k_t=cell(8,1);
 y_x=cell(8,1);

 %1 classic 1-0wLC 0 1 1
 %2 wLC renormed eco 0 1 1
 %3 wlc renormed eco alt 0 1 1
 %4 wlc renormed full 1 1 2
 %5 ML normalen(0) 0 1 1
 %6 Ml normalen 3 step(1) eco 0 1 1
 %7 ML newton alt eco(1) 0 1 1
 %8 Ml newton 1 1 2 (2) wlc renormed
 
[a_x{1},b_x{1},k_t{1},~,~,~]=Fit_wLC(lnM,W);
[a_x{2},b_x{2},k_t{2},y_x{2}]=explanation_renorm(a_x{1},b_x{1},k_t{1},g(1,:) ,norming);    
[a_x{3},b_x{3},k_t{3},y_x{3}]=explanation_renorm(a_x{1},b_x{1},k_t{1},g(1,:) ,1);       
[a_x{4},b_x{4},k_t{4},y_x{4}]=explanation_renorm(a_x{1},b_x{1},k_t{1},g ,0);    

[a_x{5},b_x{5},k_t{5}]=PoissonFit(D, E, a_startML,b_startML,k_startML); 

[a_x{6},b_x{6},k_t{6},y_x{6},~,~,~]=FitLC_Explanation_2step(lnM,g(1,:),norming,W);
[a_x{6},b_x{6},k_t{6},y_x{6}]=PoissonFit_explanation(D, E,g(1,:),a_x{6},b_x{6},k_t{6},y_x{6},norming);

[a_x{7},b_x{7},k_t{7},y_x{7}]=PoissonFit_explanation_newton(D, E,g(1,:), a_startML,b_startML,k_startML,y_startML(:,1),1);    
[a_x{8},b_x{8},k_t{8},y_x{8}]=PoissonFit_explanation_newton(D, E,g, a_x{4},b_x{4},k_t{4},y_x{4},norming);    

y_x{1}=zeros(x,n);
y_x{5}=zeros(x,n);

y_x{2}=[y_x{2},zeros(x,2)];
y_x{3}=[y_x{3},zeros(x,2)];
y_x{6}=[y_x{6},zeros(x,2)];
y_x{7}=[y_x{7},zeros(x,2)];

k_future_mean=cell(8,1);
k_future_paths=cell(8,1);
q_future_mean=cell(8,1);
q_future_paths=cell(n_paths,8);
for ii=1:8
    Model=arima(0,1,0);
%     if ii==4 || ii==8
%          Model=arima(1,1,2);
%     end
    EstM=estimate(Model,k_t{ii});
    res = infer(EstM,k_t{ii});
    k_future_mean{ii} = forecast(EstM,t_future,'Y0',k_t{ii});
    k_future_paths{ii} = simulate(EstM,t_future,'NumPaths',n_paths,'Y0',k_t{ii},'E0',res);


q_future_mean{ii}=zeros(x,t_future); 
        for i=1:x
            for j=1:t_future
               %q_future_mean(i,j)=exp(a_x(i)+b_x(i)*k_future_mean(j)); 
               %q_future_mean(i,j)=exp(a_x(i)+b_x(i)*k_future_mean(j)+y_x(i)*g_future(1,j));
    
               q_future_mean{ii}(i,j)=exp(a_x{ii}(i)+b_x{ii}(i)*k_future_mean{ii}(j)+y_x{ii}(i,:)*g_future(:,j));
             
             
            end
        end
q_future_mean{ii}=1-exp(-q_future_mean{ii});
 
    for k=1:n_paths    
       q_future_paths{k,ii}=zeros(x,t_future); 
        for i=1:x
            for j=1:t_future
               
               %q_future_paths{k}(i,j)=exp((a_x(i)+b_x(i)*k_future_paths(j,k)));
               %q_future_mean(i,j)=exp(a_x(i)+b_x(i)*k_future_mean(j)+y_x(i)*g_future_paths(1,j));
               q_future_paths{k,ii}(i,j)=exp(a_x{ii}(i)+b_x{ii}(i)*k_future_paths{ii}(j,k)+y_x{ii}(i,:)*g_future_paths{k}(:,j));
                
            end
        end
    q_future_paths{k,ii}=1-exp(-q_future_paths{k,ii});    
    end

end

    
q_M=1-exp(-M);
alpha=0.1;
i=0.009;
%X=[1 26 41 66 81]; %aged 0 25 40 65 80 
X=linspace(1,x,x);
T=[61]; %in 2016

e_Mean=cell(8,1);
e_UB=cell(8,1);
e_LB=cell(8,1);
anu_Mean=cell(8,1);
anu_UB=cell(8,1);
anu_LB=cell(8,1);

for ii=1:8
[e_Mean{ii},e_UB{ii},e_LB{ii}] = lifeexp_q2(X, T, q_M,x_min, x_max, q_future_mean{ii},q_future_paths(:,ii), alpha);
[anu_Mean{ii},anu_UB{ii},anu_LB{ii}] = ax_q2(X, T, q_M,x_min, x_max, q_future_mean{ii},q_future_paths(:,ii), alpha,i);
end
 
x_cut=71; %ie 80
hhhx=figure;
plot(linspace(0,x_cut-1,x_cut),e_Mean{1}(1:x_cut),'r')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{5}(1:x_cut),'c')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{3}(1:x_cut),'g')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{7}(1:x_cut),'m')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{6}(1:x_cut),'y')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{2}(1:x_cut),'k')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{4}(1:x_cut),'b')
hold on;
plot(linspace(0,x_cut-1,x_cut),e_Mean{8}(1:x_cut),'color','1,0.4,0')
legend('1-0-wLC_{classic}','ML-normalen_{classic}','wLC-renormed-alt_{eco}','ML-newton-alt_{eco}','ML-normalen(1)_{eco}','wLC-renormed_{eco}','wLC-renormed_{full}','ML-newton(2)_{full}')
% plot(linspace(0,x_cut-1,x_cut),e_UB{1}(1:x_cut),'r--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{2}(1:x_cut),'y--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{3}(1:x_cut),'g--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{4}(1:x_cut),'b--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{5}(1:x_cut),'c--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{6}(1:x_cut),'k--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{7}(1:x_cut),'m--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{8}(1:x_cut),'color','1,0.4,0','LineStyle','--')
% 
% plot(linspace(0,x_cut-1,x_cut),e_LB{1}(1:x_cut),'r--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{2}(1:x_cut),'y--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{3}(1:x_cut),'g--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{4}(1:x_cut),'b--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{5}(1:x_cut),'c--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{6}(1:x_cut),'k--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{7}(1:x_cut),'m--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{8}(1:x_cut),'color','1,0.4,0','LineStyle','--')

%set(gca, 'YScale', 'log')

title('e_{x,2016}')
set(hhhx,'Units','Inches');
pos = get(hhhx,'Position');
set(hhhx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhx,'ext','-dpdf','-r0')

x_cut1=25;
x_cut2=90; %ie 80
hhhx=figure;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{1}(x_cut1+1:x_cut2+1)),'r')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{5}(x_cut1+1:x_cut2+1)),'c')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{3}(x_cut1+1:x_cut2+1)),'g')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{7}(x_cut1+1:x_cut2+1)),'m')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{6}(x_cut1+1:x_cut2+1)),'y')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{2}(x_cut1+1:x_cut2+1)),'k')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{4}(x_cut1+1:x_cut2+1)),'b')
hold on;
plot(linspace(x_cut1,x_cut2,x_cut2-x_cut1+1),exp(anu_Mean{8}(x_cut1+1:x_cut2+1)),'color','1,0.4,0')
legend('1-0-wLC_{classic}','ML-normalen_{classic}','wLC-renormed-alt_{eco}','ML-newton-alt_{eco}','ML-normalen(1)_{eco}','wLC-renormed_{eco}','wLC-renormed_{full}','ML-newton(2)_{full}','location','best')
% plot(linspace(0,x_cut-1,x_cut),e_UB{1}(1:x_cut),'r--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{2}(1:x_cut),'y--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{3}(1:x_cut),'g--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{4}(1:x_cut),'b--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{5}(1:x_cut),'c--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{6}(1:x_cut),'k--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{7}(1:x_cut),'m--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_UB{8}(1:x_cut),'color','1,0.4,0','LineStyle','--')
% 
% plot(linspace(0,x_cut-1,x_cut),e_LB{1}(1:x_cut),'r--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{2}(1:x_cut),'y--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{3}(1:x_cut),'g--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{4}(1:x_cut),'b--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{5}(1:x_cut),'c--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{6}(1:x_cut),'k--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{7}(1:x_cut),'m--')
% hold on;
% plot(linspace(0,x_cut-1,x_cut),e_LB{8}(1:x_cut),'color','1,0.4,0','LineStyle','--')

%set(gca, 'YScale', 'exp')

title('exp(ä_{x,2016})')
set(hhhx,'Units','Inches');
pos = get(hhhx,'Position');
set(hhhx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhx,'anuxt','-dpdf','-r0')


ax_table=zeros(8,5);
ax_table_up=zeros(8,5);
ax_table_down=zeros(8,5);
for ii=1:8
    for j=1:5
        if j==1
            k=1;
        end
          if j==2
            k=26;
          end
          if j==3
            k=41;
          end
          if j==4
            k=66;
          end
          if j==5
            k=81;
          end
       ax_table(ii,j)= e_Mean{ii}(k);
       ax_table_up(ii,j)= e_UB{ii}(k);
       ax_table_down(ii,j)= e_LB{ii}(k);
    end
end

ax_table_sw=zeros(8,5);
ax_table_up_sw=zeros(8,5);
ax_table_down_sw=zeros(8,5);

% 1 5 3 7 6 2 4 8
ax_table_sw(1,:)=ax_table(1,:);
ax_table_sw(2,:)=ax_table(5,:);
ax_table_sw(3,:)=ax_table(3,:);
ax_table_sw(4,:)=ax_table(7,:);
ax_table_sw(5,:)=ax_table(6,:);
ax_table_sw(6,:)=ax_table(2,:);
ax_table_sw(7,:)=ax_table(4,:);
ax_table_sw(8,:)=ax_table(8,:);

ax_table_up_sw(1,:)=ax_table_up(1,:);
ax_table_up_sw(2,:)=ax_table_up(5,:);
ax_table_up_sw(3,:)=ax_table_up(3,:);
ax_table_up_sw(4,:)=ax_table_up(7,:);
ax_table_up_sw(5,:)=ax_table_up(6,:);
ax_table_up_sw(6,:)=ax_table_up(2,:);
ax_table_up_sw(7,:)=ax_table_up(4,:);
ax_table_up_sw(8,:)=ax_table_up(8,:);

ax_table_down_sw(1,:)=ax_table_down(1,:);
ax_table_down_sw(2,:)=ax_table_down(5,:);
ax_table_down_sw(3,:)=ax_table_down(3,:);
ax_table_down_sw(4,:)=ax_table_down(7,:);
ax_table_down_sw(5,:)=ax_table_down(6,:);
ax_table_down_sw(6,:)=ax_table_down(2,:);
ax_table_down_sw(7,:)=ax_table_down(4,:);
ax_table_down_sw(8,:)=ax_table_down(8,:);

fprintf([repmat('%.5f\t', 1, size(ax_table_up_sw, 2)) '\n'], ax_table_up_sw')

