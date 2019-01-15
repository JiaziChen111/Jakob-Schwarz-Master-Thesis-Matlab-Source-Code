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

BuildGDPvector;
BuildDrugData2;
alk=alk_2;

 Md_g=arima(1,1,2);   
    EstMd_g = estimate(Md_g,alk');
    h=forecast(EstMd_g,t_future,'Y0', alk');
   
%alk_2  from 1950-2005 1:56
%rauch  from 1950-2003 1:54
%g_full from 1950-2016 1:67
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
%1956-1996 (1:46)
g=[g_full(7:47)';rauch(7:47)';alk(7:47)];
cov(rauch(7:47),alk(7:47))

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

n=size(g);
n=n(1);

a_start=zeros(x,1);
b_start=ones(x,1);
k_start=ones(t,1);
y_start=zeros(x,n);

a_startML=zeros(x,1);
b_startML=ones(x,1);
k_startML=zeros(t,1);
y_startML=zeros(x,n);

a_x=cell(12,1);
b_x=cell(12,1);
k_t=cell(12,1);
y_x=cell(12,1);

%1 wlcre
%2 3steo
%3 wlc newton
%4 wlc newton D
%5 wlc newton 3step
%6 Ml newton
%7 ML normalen
%8 Ml normalen 3step
%9 Ml nnewton 3step

[a_wLC,b_wLC,k_wLC,fval_wLC,exitflag_wLC,steps_wLC]=Fit_wLC(lnM,W);
[a_x{1},b_x{1},k_t{1},y_x{1}]=explanation_renorm(a_wLC,b_wLC,k_wLC,g ,norming);

[a_x{2},b_x{2},k_t{2},y_x{2},fval_2step,exitflag_2step,steps_2steps]=FitLC_Explanation_2step(lnM,g,norming,W);
%[a_2stepD,b_2stepD,k_2stepD,y_2stepD,fval_2stepD,exitflag_2stepD,steps_2stepsD]=FitLC_Explanation_2step(lnM,g,norming,D);


[a_x{3},b_x{3},k_t{3},y_x{3},fval_ex_newton,exitflag_ex_newton]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,W,norming);

[a_x{4},b_x{4},k_t{4},y_x{4},fval_D,exitflag_D]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,D,norming);

 
[a_x{6},b_x{6},k_t{6},y_x{6},fval_ML_newton,exitflag_ML_newton]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,norming);

        
[a_x{7},b_x{7},k_t{7},y_x{7},fval_ML]=PoissonFit_explanation(D, E,g, a_startML,b_startML,k_startML,y_startML,norming);


%tic
%[a_normalen,b_normalen,k_normalen,y_normalen,fval_normalen,exitflag_normalen,steps_normalen]=FitLC_GDP_normalen(lnM,g_t,a_start,b_start,k_start,y_start,W);
%timeLC_normalen=toc;

%2steo solution
a_start=a_x{2};
b_start=b_x{2};
k_start=k_t{2};
y_start=y_x{2};

a_startML=a_x{2};
b_startML=b_x{2};
k_startML=k_t{2};
y_startML=y_x{2};


[a_x{9},b_x{9},k_t{9},y_x{9},fval_ML_newton2,exitflag_ML_newton2]=PoissonFit_explanation_newton(D, E,g, a_startML,b_startML,k_startML,y_startML,norming);
    
[a_x{8},b_x{8},k_t{8},y_x{8},fval_ML2]=PoissonFit_explanation(D, E,g, a_startML,b_startML,k_startML,y_startML,norming);

[a_x{5},b_x{5},k_t{5},y_x{5},fval_ex_newton2,exitflag_ex_newton2]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,W,norming);

a_start=a_x{1};
b_start=b_x{1};
k_start=k_t{1};
y_start=y_x{1};

a_startML=a_x{1};
b_startML=b_x{1};
k_startML=k_t{1};
y_startML=y_x{1};


[a_x{12},b_x{12},k_t{12},y_x{12},fval_ML_newton3,exitflag_ML_newton3]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,norming);
    
[a_x{11},b_x{11},k_t{11},y_x{11},fval_ML3]=PoissonFit_explanation(D, E,g, a_startML,b_startML,k_startML,y_startML,norming);

[a_x{10},b_x{10},k_t{10},y_x{10},fval_ex_newton3,exitflag_ex_newton3]=FitLC_Explanation_newton(lnM,g,a_start,b_start,k_start,y_start,W,norming);

%%%%%%%%%%%%%%%
z=12; %anzahl fttings
zz=9; % 3x3 anzahl arima
lnM_fit=cell(z,1);
for ii=1:z
    lnM_fit{ii}=zeros(x,t); 
    for i=1:x
        for j=1:t
           lnM_fit{ii}(i,j)=(a_x{ii}(i)+b_x{ii}(i)*k_t{ii}(j)+y_x{ii}(i,:)*g(:,j));
        end
    end
end
M_fit=cell(z,1);
for ii=1:z    
          M_fit{ii}=exp(lnM_fit{ii});
end

t_future=24;
n_paths=10;

k_future=cell(z,zz);

%%%%%%%%%%%%%%
 %{
 for ii=1:9   
     for ll=1:5
        if ll==1
            h=diff(k_t{ii})';
            Dt=[h(2:end);diff(diff(g))]';
            Md=varm(2,2);
            EstMd = estimate(Md,Dt);
            Dfuture=forecast(EstMd,t_future,Dt);
      
            k_future{ii,ll}=zeros(t_future,1);
            k_future{ii,ll}(1)=Dfuture(1,1)+k_t{ii}(end);        
            for i=1:t_future-1
               k_future{ii,ll}(i+1)=Dfuture(i+1,1)+k_future{ii,ll}(i,1);
            end
            %ddg=gt-2gt-1-gt-2
            %%gt=ddg+2gt-1-gt-2
            %%%need to be changed if multi dimensional
            g_future{ll,1}=zeros(t_future,1);
            g_future{ll,1}(1)=Dfuture(1,2)+2*g(end)-g(end-1);  
            g_future{ll,1}(2)=Dfuture(2,2)+2*g_future{ll,1}(1)-g(end);
            for i=2:t_future-1
                g_future{ll,1}(i+1)=Dfuture(i+1,2)+2* g_future{ll,1}(i,1)-g_future{ll,1}(i-1,1);
            end     
            
        elseif ll==2
            h=diff(k_t{ii})';
            Dt=[h(2:end);diff(diff(g))]';
            Md=varm(2,1);
            EstMd = estimate(Md,Dt);
            Dfuture=forecast(EstMd,t_future,Dt);
      
            k_future{ii,ll}=zeros(t_future,1);
            k_future{ii,ll}(1)=Dfuture(1,1)+k_t{ii}(end);        
            for i=1:t_future-1
               k_future{ii,ll}(i+1)=Dfuture(i+1,1)+k_future{ii,ll}(i,1);
            end
            %ddg=gt-2gt-1-gt-2
            %%gt=ddg+2gt-1-gt-2
            %%%need to be changed if multi dimensional
            g_future{ll,1}=zeros(t_future,1);
            g_future{ll,1}(1)=Dfuture(1,2)+2*g(end)-g(end-1);  
            g_future{ll,1}(2)=Dfuture(2,2)+2*g_future{ll,1}(1)-g(end);
            for i=2:t_future-1
                g_future{ll,1}(i+1)=Dfuture(i+1,2)+2* g_future{ll,1}(i,1)-g_future{ll,1}(i-1,1);
            end
           
       
        elseif ll==3
            Md_k=arima(0,1,0);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,ll}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
          
        elseif ll==4
            Md_k=arima(0,1,1);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,ll}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
           
        else  %ll==5
            Md_k=arima(1,1,2);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,ll}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
        end    
     end
 end   
  %%seperate arima model for g for ll=3:5 %%%%this needs to be changed if g
  %%it multiple dimensions aswel to varma or speerated arima or something

    
    Md_g=arima(4,2,2);
     EstMd_g = estimate(Md_g,g_full(1:47));
     h=forecast(EstMd_g,t_future,'Y0', g_full(1:47));
     g_future{3}=h;
     
     Md_g=arima(4,2,4);
     EstMd_g = estimate(Md_g,g_full(1:47));
     h=forecast(EstMd_g,t_future,'Y0', g_full(1:47));
     g_future{4}=h;
         
     Md_g=arima(4,2,6);
     EstMd_g = estimate(Md_g,g_full(1:47));
     h=forecast(EstMd_g,t_future,'Y0', g_full(1:47));
     g_future{5}=h;
     %}
     %%%%%%%%%%%%%%%%%%%% r
     %{
  for ii=1:9   
     for ll=1:5
        if ll==1
            Dt=[diff(k_t{ii})';diff(r_full(6:46))']';
            Md=varm(2,2);
            EstMd = estimate(Md,Dt);
            Dfuture=forecast(EstMd,t_future,Dt);
      
            k_future{ii,ll}=zeros(t_future,1);
            k_future{ii,ll}(1)=Dfuture(1,1)+k_t{ii}(end);        
            for i=1:t_future-1
               k_future{ii,ll}(i+1)=Dfuture(i+1,1)+k_future{ii,ll}(i,1);
            end
            r_future=zeros(t_future,1);
            r_future(1)=Dfuture(1,2)+r_full(46);        
            for i=1:t_future-1
               r_future(i+1)=Dfuture(i+1,2)+r_future(i,1);
            end
            %ddg=gt-2gt-1-gt-2
            %%gt=ddg+2gt-1-gt-2
            %%%need to be changed if multi dimensional
            g_future{ll,1}=zeros(t_future,1);
            g_future{ll,1}(1)=(1+r_future(1))*exp(g_full(47));  
            
            for i=1:t_future-1
                g_future{ll,1}(i+1)=(1+r_future(i))*g_future{ll,1}(i);
            end     
             g_future{ll,1}=log( g_future{ll,1});
        elseif ll==2
            Dt=[diff(k_t{ii})';diff(r_full(6:46))']';
            Md=varm(2,1);
            EstMd = estimate(Md,Dt);
            Dfuture=forecast(EstMd,t_future,Dt);
      
            k_future{ii,ll}=zeros(t_future,1);
            k_future{ii,ll}(1)=Dfuture(1,1)+k_t{ii}(end);        
            for i=1:t_future-1
               k_future{ii,ll}(i+1)=Dfuture(i+1,1)+k_future{ii,ll}(i,1);
            end
            %ddg=gt-2gt-1-gt-2
            %%gt=ddg+2gt-1-gt-2
            %%%need to be changed if multi dimensional
            r_future=zeros(t_future,1);
            r_future(1)=Dfuture(1,2)+r_full(46);        
            for i=1:t_future-1
               r_future(i+1)=Dfuture(i+1,2)+r_future(i,1);
            end
            %ddg=gt-2gt-1-gt-2
            %%gt=ddg+2gt-1-gt-2
            %%%need to be changed if multi dimensional
               g_future{ll,1}=zeros(t_future,1);
            g_future{ll,1}(1)=(1+r_future(1))*exp(g_full(47));  
            
            for i=1:t_future-1
                g_future{ll,1}(i+1)=(1+r_future(i))*g_future{ll,1}(i);
            end     
             g_future{ll,1}=log( g_future{ll,1});
           
       
        elseif ll==3
            Md_k=arima(0,1,0);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,ll}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
          
        elseif ll==4
            Md_k=arima(0,1,1);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,ll}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
           
        else  %ll==5
            Md_k=arima(1,1,2);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,ll}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
        end    
     end
 end  
     
     Md_g=arima(1,1,3);
     EstMd_g = estimate(Md_g,r_full(1:46));
     h=forecast(EstMd_g,t_future,'Y0', r_full(1:46));
     
     g_future{3,1}=zeros(t_future,1);
     g_future{3}(1)=(1+h(1))*exp(g_full(47));
     for i=2:length(g_future{3})
         g_future{3}(i)=(1+h(i))*g_future{3}(i-1);
     end
      g_future{3}= log(g_future{3});
     
     Md_g=arima(2,1,2);
     EstMd_g = estimate(Md_g,r_full(1:46));
     h=forecast(EstMd_g,t_future,'Y0', r_full(1:46));
      g_future{4,1}=zeros(t_future,1);
     g_future{4}(1)=(1+h(1))*exp(g_full(47));
     for i=2:length(g_future{4})
         g_future{4}(i)=(1+h(i))*g_future{4}(i-1);
     end
            g_future{4}= log(g_future{4}); 
            
     Md_g=arima(1,1,2);
     EstMd_g = estimate(Md_g,r_full(1:46));
     h=forecast(EstMd_g,t_future,'Y0', r_full(1:46));
      g_future{5,1}=zeros(t_future,1);
     g_future{5}(1)=(1+h(1))*exp(g_full(47));
     for i=2:length(g_future{5})
         g_future{5}(i)=(1+h(i))*g_future{5}(i-1);
     end
         g_future{5}= log(g_future{5});
%%%%%%%%%%%%%%%%%
    
   gforecast=figure;
             plot(linspace(1950,1950+67,67),g_full(1:end),'g')
             hold on
             plot(linspace(42+1955,65+1955,24),g_future{1}(1:24),'r')
             plot(linspace(42+1955,65+1955,24),g_future{2}(1:24),'m')
             plot(linspace(42+1955,65+1955,24),g_future{3}(1:24),'b')
             plot(linspace(42+1955,65+1955,24),g_future{4}(1:24),'c')
             plot(linspace(42+1955,65+1955,24),g_future{5}(1:24),'k')
             title('forcasting g_t')
             legend('historical','varma(2,0)','varma(1,0)','arima(1,1,3)','arima(2,1,2)','arima(1,1,2)' ,'location','best');
            set(gforecast,'Units','Inches');
pos = get(gforecast,'Position');
set(gforecast,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gforecast,'gforecast3','-dpdf','-r0')

             %varma 4,0 varma 0 0 and arima 0 2 0 almost idetical, varma 3 0 even. 
%varma 1 0 slighly better than varma 2 0
%stronger dorwards slope
 %}
  for ii=1:z   
     k=1;
        for i=0:2
            for j=0:2
            Md_k=arima(i,1,j);
            EstMd_k = estimate(Md_k,k_t{ii});
            k_future{ii,k}=forecast(EstMd_k,t_future,'Y0',k_t{ii});
            k=k+1;
            end
        end    
     
 end  
     %g_future=cell(5,1);
     g_future=zeros(n,t_future);
     Md_g=arima(1,1,2);
     EstMd_g = estimate(Md_g,r_full(1:46));
     h=forecast(EstMd_g,t_future,'Y0', r_full(1:46));
     g_future(1,1)=(1+h(1))*exp(g_full(47));
     for i=2:t_future
         g_future(1,i)=(1+h(i))*g_future(1,i-1);
     end
         g_future(1,:)= log(g_future(1,:))';
    
%      %alk
%      pBIC=0;
%      qBIC=0;
%      pAIC=0;
%      qAIC=0;
%      Md=arima(0,1,0);  
%              EstMd = estimate(Md,alk');
%              results = summarize(EstMd);
%      AIC=results.AIC;
%      BIC=results.BIC;
%      for i=0:5
%          for j=0:5
%              i
%              j
%          if (i==5&& j==4) ||  (i==5&& j==5)  
%              continue
%          end
%              Md=arima(i,1,j);  
%              EstMd = estimate(Md,alk');
%              results = summarize(EstMd);
%         if results.AIC<AIC
%             AIC=results.AIC;
%             pAIC=i;
%             qAIC=j;
%         end
%         if results.BIC<BIC
%             BIC=results.BIC;
%             pBIC=i;
%             qBIC=j;
%         end
%          end
%     end
    
    %best 
    Md_g=arima(1,1,2);   
    EstMd_g = estimate(Md_g,alk(1:46)');
    h=forecast(EstMd_g,t_future,'Y0', alk(1:46)');
    g_future(2,:)=h;
    
  
%   alky=figure;
%  plot(linspace(1950,1990,41),alk_2(1:41),'k');
%  hold on;
%  plot(linspace(1991,2005,15),alk_2(42:end),'r');
%  hold on
%  plot(linspace(1997,2020,24),h,'b--');
%  hold on;
%  Md_g=arima(1,1,1);   
%     EstMd_g = estimate(Md_g,alk(1:46)');
%     h=forecast(EstMd_g,t_future,'Y0', alk(1:46)');
%     g_future(2,:)=h;
%  plot(linspace(1997,2020,24),h,'c--');
%  title('Daily alkohol consumption: g^{(3)}_t');
%  legend('West Germany','Germany','arima(1,1,2)','arima(1,1,1)')
%  set(alky,'Units','Inches');
%  pos = get(alky,'Position');
%  set(alky,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%  print(alky,'alky','-dpdf','-r0')

    
    
    %rauch
    %1950 1955 165 1975 1985 1996
    h=interp1(linspace(1950,1996,54),rauch',linspace(1997,2020,24),'linear','extrap');
    %plot(x_org,rauch_m_org,'o',linspace(1997,2020,24),h,':.');
    g_future(3,:)=h;
    
%     uu=figure;
%     subplot(3,1,1)
%     h=interp1(linspace(1950,1996,54),rauch_m',linspace(1997,2020,24),'linear','extrap');
%     plot(x_org,rauch_m_org,'o',linspace(1997,2020,24),h,':.');
%     title ('males smoking percentage')
%     subplot(3,1,2)
%     h=interp1(linspace(1950,1996,54),rauch_f',linspace(1997,2020,24),'linear','extrap');
%     plot(x_org,rauch_f_org,'o',linspace(1997,2020,24),h,':.');
%     title ('females smoking percentage')
%     subplot(3,1,3)
%     h=interp1(linspace(1950,1996,54),rauch_ges',linspace(1997,2020,24),'linear','extrap');
%     plot(x_org,rauch_ges_org,'o',linspace(1997,2020,24),h,':.');
%     title ('total smoking percentage')
%     set(uu,'Units','Inches');
%     pos = get(uu,'Position');
%     set(uu,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(uu,'uu','-dpdf','-r0')

    
    
lnM_future=cell(z,zz);
for ll=1:zz
    for ii=1:z
        lnM_future{ii,ll}=zeros(x,t_future); 
        for i=1:x
            for j=1:t_future
               lnM_future{ii,ll}(i,j)=(a_x{ii}(i)+b_x{ii}(i)*k_future{ii,ll}(j)+y_x{ii}(i,:)*g_future(:,j));
            end
        end
    end
end
M_future=cell(z,zz);
for ii=1:z   
    for ll=1:zz
          M_future{ii,ll}=exp(lnM_future{ii,ll});
    end
end

ErrM=zeros(z,zz);
for ll=1:zz
    for ii=1:z
        ErrM(ii,ll)= errorfkt3( M_full(1:x,42:end),M_future{ii,ll}(:,1:19));
    end
end
%ErrM=ErrM([1,7,8,11,12],:);
%ErrM=ErrM';
ErrM=ErrM*10^5;

fprintf([repmat('%.5f\t', 1, size(ErrM, 2)) '\n'], ErrM')

%%%%% classic
 Md_k=arima(0,1,1);
 EstMd_k = estimate(Md_k,k_wLC);
            k_wLC_future=forecast(EstMd_k,t_future,'Y0',k_wLC);

M_wLC_future=zeros(x,t_future); 
        for i=1:x
            for j=1:t_future
               M_wLC_future(i,j)=exp((a_wLC(i)+b_wLC(i)*k_wLC_future(j)));
            end
        end
index42=M_full(1:x,42:end)==42;

Merr=cell(z,zz);
for ll=1:zz
    for ii=1:z
        Merr{ii,ll}=zeros(x,19);
        Merr{ii,ll}= (M_future{ii,ll}(:,1:19)./M_full(1:x,42:end))-1;
        Merr{ii,ll}(index42)=0;
    end
end

%gdp
Merrclassic= (M_wLC_future(:,1:19)./M_full(1:x,42:end))-1 ;
Merrclassic(index42)=0;

g_t=g(1,:);
[a_start,b_start,k_start,y_start,~,~,~]=FitLC_Explanation_2step(lnM,g_t,norming,W);
[a_gdpML,b_gdpML,k_gdpML,y_gdpML,~]=PoissonFit_explanation_newton(D, E,g_t, a_start,b_start,k_start,y_start,norming);
[a_gdp,b_gdp,k_gdp,y_gdp]=explanation_renorm(a_wLC,b_wLC,k_wLC,g_t ,norming);

Md_k=arima(0,1,1);
            EstMd_k = estimate(Md_k,k_gdpML);
            k_gdpML_future=forecast(EstMd_k,t_future,'Y0',k_gdpML);
Md_k=arima(1,1,2);
            EstMd_k = estimate(Md_k,k_gdp);
            k_gdp_future=forecast(EstMd_k,t_future,'Y0',k_gdp);

 M_gdp_future=zeros(x,t_future); 
 M_gdpML_future=zeros(x,t_future);
        for i=1:x
            for j=1:t_future
               M_gdpML_future(i,j)=exp(a_gdpML(i)+b_gdpML(i)*k_gdpML_future(j)+y_gdpML(i)*g_future(1,j));
               M_gdp_future(i,j)=exp(a_gdp(i)+b_gdp(i)*k_gdp_future(j)+y_gdp(i)*g_future(1,j));
            end
        end
index42=M_full(1:x,42:end)==42;

Merr_gdp=(M_gdp_future(:,1:19)./M_full(1:x,42:end))-1;
Merr_gdp(index42)=0;
  
Merr_gdpML=(M_gdpML_future(:,1:19)./M_full(1:x,42:end))-1;
Merr_gdpML(index42)=0;
      
        
% errq=figure;
% 
% subplot(2,2,1)
% xscale = [1997 2015];
% yscale = [0 110];
% imagesc(xscale,yscale,Merr{1,4})
% colorbar
% caxis([-1 1]);
% xlabel('year t')
% ylabel('age x')
% title('wLC_{renormed}');
% 
% subplot(2,2,2)
% imagesc(xscale,yscale,Merrclassic)
% colorbar
% caxis([-1 1]);
% xlabel('year t')
% ylabel('age x')
% title('wLC_{classic}');
% 
% subplot(2,2,3)
% xscale = [1997 2015];
% yscale = [0 110];
% imagesc(xscale,yscale,Merr{6,4})
% colorbar
% caxis([-1 1]);
% xlabel('year t')
% ylabel('age x')
% title('ML_{newton}');
% 
% subplot(2,2,4)
% imagesc(xscale,yscale,Merr{2,2})
% colorbar
% caxis([-1 1]);
% xlabel('year t')
% ylabel('age x')
% title('3step');
% 
% suptitle('relative forecast error colormap')
% 
% set(errq,'Units','Inches');
% pos = get(errq,'Position');
% set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(errq,'errgdpforecast2','-dpdf','-r0')

%%%%%%%%%%norming ML solution to compare
% a_start=zeros(x,1);
% b_start=ones(x,1);
% k_start=ones(t,1);
% y_start=zeros(x,n);

%a_xN=cell(9,1);
%b_xN=cell(9,1);
%k_tN=cell(9,1);
%y_xN=cell(9,1);
%[a_xN{6},b_xN{6},k_tN{6},y_xN{6},fval_ML_newtonN,exitflag_ML_newtonN]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,norming);
%[a_xN,b_xN,k_tN,y_xN,fval_ML_newtonN,exitflag_ML_newtonN]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,1);
%Md_k=arima(1,1,1);
%       EstMd_k = estimate(Md_k,k_tN);
%            k_futureN=forecast(EstMd_k,t_future,'Y0',k_tN);

% [a_xN,b_xN,k_tN,y_xN,fval_ML_newtonN,exitflag_ML_newtonN]=PoissonFit_explanation_newton(D, E,g, a_start,b_start,k_start,y_start,1);
% [a_renormN,b_renormN,k_renormN,y_renormN]=explanation_renorm(a_wLC,b_wLC,k_wLC,g ,1);
% Md_k=arima(0,1,1);
%        EstMd_k = estimate(Md_k,k_tN);
%             k_futureN=forecast(EstMd_k,t_future,'Y0',k_tN);
%        EstMd_k = estimate(Md_k,k_tN);
%             k_future_renormN=forecast(EstMd_k,t_future,'Y0',k_renormN);
% 
% M_futureN=zeros(x,t_future);  
%         for i=1:x
%             for j=1:t_future
%                M_futureN(i,j)=exp((a_xN(i)+b_xN(i)*k_futureN(j)+y_xN(i,:)*g_future{1}(:,j)));
%             end
%         end
% M_future_renormN=zeros(x,t_future);  
%         for i=1:x
%             for j=1:t_future
%                M_future_renormN(i,j)=exp((a_renormN(i)+b_renormN(i)*k_future_renormN(j)+y_renormN(i,:)*g_future{1}(:,j)));
%             end
%         end
%            
            
            
% 
% index42=M_full(1:x,42:end)==42;
% %wlc renormed 212
% MerrrMLN= (M_futureN(:,1:19)./M_full(1:x,42:end))-1;
% MerrrenormN= (M_future_renormN(:,1:19)./M_full(1:x,42:end))-1 ;
% 
% MerrenormN(index42)=0;
% MerrMLN(index42)=0;          

errq=figure;
subplot(2,2,1)
xscale = [1997 2015];
yscale = [0 100];
imagesc(xscale,yscale,Merr{1,6})
colorbar
caxis([-1 1]);
xlabel('year t')
ylabel('age x')
title('wLC-renormed_{full}');

subplot(2,2,2)
xscale = [1997 2015];
yscale = [0 100];
imagesc(xscale,yscale,Merrclassic)
colorbar
caxis([-1 1]);
xlabel('year t')
ylabel('age x')
title('1-0-wLC_{classic}');

subplot(2,2,3)
xscale = [1997 2015];
yscale = [0 100];
imagesc(xscale,yscale,Merr_gdp)
colorbar
caxis([-1 1]);
xlabel('year t')
ylabel('age x')
title('wLC-renormed_{eco}');

subplot(2,2,4)
imagesc(xscale,yscale,Merr_gdpML)
colorbar
caxis([-1 1]);
xlabel('year t')
ylabel('age x')
title('ML-newton(1)_{eco}');

suptitle('relative forecast error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errgdpforecast3rauch','-dpdf','-r0')

%%%%%%%%%%
x_cut=x   ;   
hhhx=figure;

subplot(2,2,1);
plot(linspace(1-1,x_cut-1 ,x_cut ),a_x{1}(1:x_cut ),'r');
hold on;
% plot(linspace(1,x_cut-1 ,x_cut ),a_wLC(1:x_cut ),'b');
% hold on;
plot(linspace(1,x_cut-1 ,x_cut ),a_gdpML(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut-1 ,x_cut ),a_gdp(1:x_cut ),'k');

title('a_x');

subplot(2,2,2);
plot(linspace(1-1,x_cut-1 ,x_cut ),b_x{1}(1:x_cut ),'r');
hold on;
% plot(linspace(1,x_cut-1 ,x_cut ),b_wLC(1:x_cut ),'b');
% hold on;
plot(linspace(1,x_cut-1 ,x_cut ),b_gdpML(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut-1 ,x_cut ),b_gdp(1:x_cut ),'k');
%legend({'wLC-renormed_{full}','1-0-wLC_{classic}','ML-newton(1)_{eco}','wLC-renormed_{eco}'},'location','best');
legend({'wLC-renormed_{full}','ML-newton(1)_{eco}','wLC-renormed_{eco}'},'location','best');

title('b_x');

subplot(2,2,3);
plot(linspace(1,t ,t ),k_t{1},'r');
hold on;
% plot(linspace(1,t ,t ),k_wLC,'b');
% hold on;
plot(linspace(1,t ,t ),k_gdpML,'b');
hold on;
plot(linspace(1,t ,t ),k_gdp,'k');

plot(linspace(t+1,t+t_future ,t_future ),k_future{1,6},'r--');
hold on;
% plot(linspace(t+1,t+t_future ,t_future),k_wLC_future,'b--');
% hold on;
plot(linspace(t+1,t+t_future ,t_future),k_gdpML_future,'b--');
hold on;
plot(linspace(t+1,t+t_future ,t_future),k_gdp_future,'k--');


title('k_t');

subplot(2,2,4);
plot(linspace(1-1,x_cut-1 ,x_cut ),y_x{1}(1:x_cut,1 ),'r');
hold on;
plot(linspace(1,x_cut-1 ,x_cut ),y_gdpML(1:x_cut ),'b');
hold on;
plot(linspace(1,x_cut-1 ,x_cut ),y_gdp(1:x_cut ),'k');
title('y^{(1)}_x');
set(hhhx,'Units','Inches');
pos = get(hhhx,'Position');
set(hhhx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhx,'axbxktyxforecast100rauch','-dpdf','-r0')

hhhx=figure;

subplot(2,1,1);
plot(linspace(1-1,x_cut-1 ,x_cut ),y_x{1}(1:x_cut,2 ),'r');

title('y^{(2)}_x');



subplot(2,1,2);
plot(linspace(1-1,x_cut-1 ,x_cut ),y_x{1}(1:x_cut,3 ),'r');

title('y^{(3)}_x');

set(hhhx,'Units','Inches');
pos = get(hhhx,'Position');
set(hhhx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhx,'axbxktyxforecast100rauch2','-dpdf','-r0')

%%%%%%%%%% jenkins
hhhx=figure;
subplot(2,2,1);
plot(k_t{1})
title('k_t');
subplot(2,2,2);
autocorr(k_t{1})
title(' ');
subplot(2,2,3);
plot(diff(k_t{1}))
title('dk_t');
subplot(2,2,4);
autocorr(diff(k_t{1}))
title(' ');

% title('y^{(3)}_x');

set(hhhx,'Units','Inches');
pos = get(hhhx,'Position');
set(hhhx,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhhx,'jenkins','-dpdf','-r0')
