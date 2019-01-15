clear all;
%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=110; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 40; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42
%load data with above parameters
BuildCentralDeathMatrix; 
%%%Cairns Model%%%
    x_max_C=x_max; %for a new maximum age group, i .e. 100 instead of 110
    x_min_C= x_min; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1, A_2,mu,C]=q_RegressionFit(q,x_min_C,x_max_C);

%set some data parameters
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =3; % chose sex 3 4 5 for female male total
x_max=110; %new maximum age group 100+, set to 104 to cut of all -1,0 in M
x_min= 70; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included for LC like models
t_max=2015; %last year of data included for LC like models
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
replace42=1; %if 1 then all values <= in M are replaced with 42
%load data with above parameters
BuildCentralDeathMatrix; 
%%%Cairns Model%%%
    x_max_C=x_max; %for a new maximum age group, i .e. 100 instead of 110
    x_min_C= x_min; %for a new minimum included age group, i .e. 60 instead of 0
    t_min_C=t_min; %for first year of data to include. i.e. 1960 instead of 1956
    t_max_C=t_max; %for last year of data to include. i.e. 2010 instead of 2015
    BuildDeathProbMatrix;% creates Matrix q for observed Deathprobabilities
    %changes needed to adapt for other data inside aboves script:
       %fltper_1x1.txt
       %mltper_1x1.txt
       %bltper_1x1.txt
       %note: sex is used from BuildCentralDeathMatrix
    [A_1_70, A_2_70,mu_70,C_70]=q_RegressionFit(q,x_min_C,x_max_C);



q_Plots;

h1=figure;
plot(A_1,A_2);
hold on
plot(A_1_70,A_2_70,'r');
xlabel('A_1(t)');
ylabel('A_2(t)');
title('A(t))');
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1,'At','-dpdf','-r0')

%%%%%%%
h2=figure;
dA=[diff(A_1),diff(A_2)];
%plot(dA(:,1),dA(:,2),'*');
hold on;
% Calculate the eigenvectors and eigenvalues
covariance = C*C';
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
avg = mu;

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
print(h2,'dAt','-dpdf','-r0')
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
print(hA1,'A1','-dpdf','-r0')

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
print(hA2,'A2','-dpdf','-r0')

%%%%%%%%%%
errq=figure;
subplot(2,1,1)
Merr= q- q_C ;

xscale = [t_min_C t_max_C];
yscale = [x_min_C x_max_C];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.01 0.01]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,1,2)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.001 0.001]);
xlabel('year t')
ylabel('age x')
title('')

suptitle('Fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errq70','-dpdf','-r0')

%%%%%%%%%%
index42=M(:,:)==42;
%qx=1-exp(-mx),mx=-log(1-qx),
M_Cairns=-log(1-q_C);
errq=figure;
subplot(2,2,1)
Merr= M- M_Cairns ;
Merr(index42)=0;
xscale = [t_min_C t_max_C];
yscale = [x_min_C x_max_C];
imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.5 1]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,2)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.05 0.075]);
xlabel('year t')
ylabel('age x')
title('')

subplot(2,2,3)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.005 0.0075]);
xlabel('year t')
ylabel('age x')

subplot(2,2,4)

imagesc(xscale,yscale,Merr)
colorbar
caxis([-0.0005 0.00075]);
xlabel('year t')
ylabel('age x')
suptitle('Fitting error colormap')

set(errq,'Units','Inches');
pos = get(errq,'Position');
set(errq,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(errq,'errq270','-dpdf','-r0')
