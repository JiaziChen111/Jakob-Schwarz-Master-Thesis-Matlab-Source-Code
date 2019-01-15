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
q_M=1-exp(-M_full);
%%%Cairns Model%%%
    x_max_C=110; %for a new maximum age group, i .e. 100 instead of 110
    x_min_C= 40; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
   q=zeros(x_C,t_C);
    q(1:x_C,1:t_C)=q_M(x_min_C+1:x_max_C+1,(t_min_C-t_start_C+1):t_full_C-(t_end_C-t_max_C));

    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
       
 [A_1, A_2]=q_RegressionFit2(q,x_min_C,x_max_C);

        x_max_C_70=110; %for a new maximum age group, i .e. 100 instead of 110
        x_min_C_70= 70; %for a new minimum included age group, i .e. 60 instead of 0
        x_C_70=x_max_C_70-x_min_C_70+1; %number of age groups used
        [A_1_70, A_2_70]=q_RegressionFit2(q(31:end,:),x_min_C_70,x_max_C_70);

%%%%%%%%%%%%%%%%%%%%
q_C=zeros(x_C,t_C);
    for i=1:x_C
        for j=1:t_C
            h=A_1(j)+A_2(j)*((x_min_C+i-1)+(j-1));
            q_C(i,j)=exp(h)/(1+exp(h));
        end
    end
q_C_70=zeros(x_C_70,t_C);
    for i=1:x_C_70
        for j=1:t_C
            h=A_1_70(j)+A_2_70(j)*((x_min_C_70+i-1)+(j-1));
            q_C_70(i,j)=exp(h)/(1+exp(h));
        end
    end

 [a_wLC,b_wLC,k_wLC,steps_wLC,fval_wLC]=Fit_wLC(lnM);
  
%%%%%%%
M_wLC=zeros(x,t);
    for i=1:x
        for j=1:t
            M_wLC(i,j)=exp(a_wLC(i)+b_wLC(i)*k_wLC(j));
        end
    end  
q_wLC=1-exp(-M_wLC);
    
%%%%%%%%%%%%%  
qC=figure('name','q(x,t) Fit with Cairns Model');
ax1 = subplot(4,1,1);
%plots development over time for a specific age
age=40; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);

plot(time,q_M(age+1,:),'g')
hold on
plot(time,q_C(age-x_min_C+1,:),'b')
hold on
plot(time,q_wLC(age+1,:),'c')
h=num2str(age);

title(ax1,['q(x,t) over time for fixed age ' h])


ax9 = subplot(4,1,2);
%plots development over time for a specific age
age=70; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);

plot(time,q_M(age+1,:),'g')
hold on
plot(time,q_C(age-x_min_C+1,:),'b')
hold on
plot(time,q_C_70(age-x_min_C_70+1,:),'r')
hold on
plot(time,q_wLC(age+1,:),'c')
h=num2str(age);

title(ax9,['q(x,t) over time for fixed age ' h])


axt = subplot(4,1,3);
%plots development over time for a specific age
age=95; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);

plot(time,q_M(age+1,:),'g')
hold on
plot(time,q_C(age-x_min_C+1,:),'b')
hold on
plot(time,q_C_70(age-x_min_C_70+1,:),'r')
hold on
plot(time,q_wLC(age+1,:),'c')
h=num2str(age);

title(axt,['q(x,t) over time for fixed age ' h])
%plots development over time for a specific age
axy=subplot(4,1,4);
age=110; %1:x_C enstpricht x_min:x_max %%%%%%%%%%%%%%choose
time=linspace(t_min_C,t_max_C,t_C);
index42=M_full==42;
index42=logical((index42-1)*(-1));
plot(time(index42(age+1,:)),q_M(age+1,index42(age+1,:)),'g')
hold on
plot(time,q_C(age-x_min_C+1,:),'b')
hold on
plot(time,q_C_70(age-x_min_C_70+1,:),'r')
hold on
plot(time,q_wLC(age+1,:),'c')
h=num2str(age);

title(axy,['q(x,t) over time for fixed age ' h])

set(qC,'Units','Inches');
pos = get(qC,'Position');
set(qC,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(qC,'qCaM','-dpdf','-r0')

ppp=figure;
%plots development over age for a specific year
year=1; %beetween 1:t, i.e. t_min-t_max %%%%%%%%%%%%%% choose

h=num2str(1956);
h1x=num2str(2015);

ages=linspace(40,110,x_C);
ages_70=linspace(70,110,x_C_70);
index42=M_full==42;
index42=logical((index42-1)*(-1));
Q=q_M(ages+1,year);
index42=index42(41:end,:);
plot(ages(index42(:,year)),log(Q(index42(:,year))),'g')
hold on
plot(ages,log(q_C(:,year)),'b')
hold on
plot(ages_70,log(q_C_70(:,year)),'r')
hold on
plot(ages,log(q_wLC(ages+1,year)),'c')
hold on 
year=60;
Q=q_M(ages+1,year);
plot(ages(index42(:,year)),log(Q(index42(:,year))),'g')
hold on
plot(ages,log(q_C(:,year)),'b')
hold on
plot(ages_70,log(q_C_70(:,year)),'r')
hold on
plot(ages,log(q_wLC(ages+1,year)),'c')
title('ln q(x,t) over ages in year 1956  and 2015 (lower)')
legend('Historical (implied)','Cairns','Cairns70','wLC','Location','Best')
set(ppp,'Units','Inches');
pos = get(ppp,'Position');
set(ppp,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(ppp,'qCbM','-dpdf','-r0')

%%%%%%%
h2=figure;
dA=[diff(A_1),diff(A_2)];
%plot(dA(:,1),dA(:,2),'*');
hold on;
% Calculate the eigenvectors and eigenvalues
V=cov(dA);
covariance = V;
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(dA);

% Get the 95% confidence interval error ellipse
chisquare_val = sqrt(5.991);

theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'r--')
hold on;
data=dA;
plot(data(:,1), data(:,2), '.');
mindata = min(min(data));
maxdata = max(max(data));

hold on;

% Plot the eigenvectors
%quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
%quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);

% Get the 99% confidence interval error ellipse
chisquare_val = sqrt(9.210);

theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'k-.')
hold on;
data=dA;
plot(data(:,1), data(:,2), '.');
mindata = min(min(data));
maxdata = max(max(data));

hold on;
% Get the 90% confidence interval error ellipse
chisquare_val = sqrt(4.605);

theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'g')
hold on;
data=dA;
plot(data(:,1), data(:,2), '.');
mindata = min(min(data));
maxdata = max(max(data));

hold on;

title('A(t)-A(t-1)');
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2,'dAtM','-dpdf','-r0')
%%%%%%%%%

hA1=figure;
plot(linspace(t_min_C,t_max_C,t_C),A_1);
hold on;
plot(linspace(t_min_C,t_max_C,t_C),A_1_70,'r');
xlabel('time t');
ylabel('A_1');
title('A_1(t))');
set(hA1,'Units','Inches');
pos = get(hA1,'Position');
set(hA1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hA1,'A1M','-dpdf','-r0')

hA2=figure;
plot(linspace(t_min_C,t_max_C,t_C),A_2);
hold on;
plot(linspace(t_min_C,t_max_C,t_C),A_2_70,'r');
xlabel('time t');
ylabel('A_2');
title('A_2(t))');
set(hA2,'Units','Inches');
pos = get(hA2,'Position');
set(hA2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hA2,'A2M','-dpdf','-r0')

%%%%%%%%%%
index42=M==42;
index42=index42(41:end,:);
%qx=1-exp(-mx),mx=-log(1-qx),
errq=figure;
subplot(2,2,1)
Merr= q_M(41:end,:)- q_C ;
Merr(index42)=0;
xscale = [t_min_C t_max_C];
yscale = [x_min_C x_max_C];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.3 0.3]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.1 0.1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.05]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')

suptitle('Fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errq2M','-dpdf','-r0')
%%%%%%%%%%%
[a_D,b_D,k_D,steps_D,fval_D]=Fit_wLC(lnM,D);
[a_D2,b_D2,k_D2]=LeeCarterFit(D,E,a_D,b_D,k_D);

%custom weights
%w_x,t=t+x_max-x 
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
 [a_wLC2,b_wLC2,k_wLC2]=LeeCarterFit(D,E,a_wLC,b_wLC,k_wLC);
 
  [a_ML,b_ML,k_ML,fval_ML]=PoissonFit(D, E);
  [a_ML2,b_ML2,k_ML2]=LeeCarterFit(D,E,a_ML,b_ML,k_ML);
  
%%%%%%%
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
M_D=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D(i,j)=exp(a_D(i)+b_D(i)*k_D(j));
        end
    end    
    M_D2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_D2(i,j)=exp(a_D2(i)+b_D2(i)*k_D2(j));
        end
    end    
  M_W=zeros(x,t);
    for i=1:x
        for j=1:t
            M_W(i,j)=exp(a_W(i)+b_W(i)*k_W(j));
        end
    end
M_W2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_W2(i,j)=exp(a_W2(i)+b_W2(i)*k_W2(j));
        end
    end
M_ML=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML(i,j)=exp(a_ML(i)+b_ML(i)*k_ML(j));
        end
    end
M_ML2=zeros(x,t);
    for i=1:x
        for j=1:t
            M_ML2(i,j)=exp(a_ML2(i)+b_ML2(i)*k_ML2(j));
        end
    end
    
q_wLC=1-exp(-M_wLC);
q_wLC2=1-exp(-M_wLC2);
q_W=1-exp(-M_W);
q_W2=1-exp(-M_W2);
q_D=1-exp(-M_D);
q_D2=1-exp(-M_D2);
q_ML=1-exp(-M_ML);
q_ML2=1-exp(-M_ML2);

ErrM=zeros(9,5);
[ MaxErr,MinErr,MeanErr,Err2,R2,R_C]=errorfkt_q( q,q_C );
ErrM(1,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_wLC]=errorfkt_q( q,q_wLC(x_min_C+1:end,:) );
ErrM(2,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_wLC2]=errorfkt_q( q,q_wLC2(x_min_C+1:end,:) );
ErrM(3,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_D]=errorfkt_q( q,q_D(x_min_C+1:end,:) );
ErrM(4,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_D2]=errorfkt_q( q,q_D2(x_min_C+1:end,:) );
ErrM(5,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_W]=errorfkt_q( q,q_W(x_min_C+1:end,:) );
ErrM(6,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_W2]=errorfkt_q( q,q_W2(x_min_C+1:end,:) );
ErrM(7,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_ML]=errorfkt_q(q,q_ML(x_min_C+1:end,:) );
ErrM(8,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_ML2]=errorfkt_q( q,q_ML2(x_min_C+1:end,:) );
ErrM(9,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

fprintf([repmat('%.4f\t', 1, size(ErrM, 2)) '\n'], ErrM')

ErrM_70=zeros(9,5);
[ MaxErr,MinErr,MeanErr,Err2,R2,R_C]=errorfkt_q( q(31:end,:),q_C_70 );
ErrM_70(1,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_wLC]=errorfkt_q( q(31:end,:),q_wLC(x_min_C_70+1:end,:) );
ErrM_70(2,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_wLC2]=errorfkt_q( q(31:end,:),q_wLC2(x_min_C_70+1:end,:) );
ErrM_70(3,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_D]=errorfkt_q( q(31:end,:),q_D(x_min_C_70+1:end,:) );
ErrM_70(4,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_D2]=errorfkt_q( q(31:end,:),q_D2(x_min_C_70+1:end,:) );
ErrM_70(5,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_W]=errorfkt_q( q(31:end,:),q_W(x_min_C_70+1:end,:) );
ErrM_70(6,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_W2]=errorfkt_q( q(31:end,:),q_W2(x_min_C_70+1:end,:) );
ErrM_70(7,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_ML]=errorfkt_q(q(31:end,:),q_ML(x_min_C_70+1:end,:) );
ErrM_70(8,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

[ MaxErr,MinErr,MeanErr,Err2,R2,R_ML2]=errorfkt_q( q(31:end,:),q_ML2(x_min_C_70+1:end,:) );
ErrM_70(9,:)=[ MaxErr,MinErr,MeanErr,Err2,R2];

fprintf([repmat('%.4f\t', 1, size(ErrM_70, 2)) '\n'], ErrM_70')

xzf=figure('name','FVU_x');
 [ MaxErr,MinErr,MeanErr,Err2,R2,R_C_70]=errorfkt_q( q(31:end,:),q_C_70 );
  
   plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_C,'b')
   hold on
   plot(linspace(70,110,x_C_70),R_C_70,'r')
   hold on
   plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_wLC,'c')
   hold on
    plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_D,'k')
   hold on
  plot(linspace(x_min_C,x_min_C+x_C-1,x_C),R_ML,'y')
 
    ylabel('FVU_x')
    xlabel('x')
   
  legend({'Cairns','Cairns70','1-0-weights','Death-weights,','ML'},'Location','Best')  
     set(xzf,'Units','Inches');
pos = get(xzf,'Position');
set(xzf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(xzf,'FVU_q2M','-dpdf','-r0')
     
    