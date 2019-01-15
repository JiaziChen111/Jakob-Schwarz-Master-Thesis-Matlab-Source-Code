%changes to adapt for other data see Master to set:

%{
DeathRateData= importdata('Mx_1x1.txt');
E= importdata('Exposures_1x1.txt');
D= importdata('Deaths_1x1.txt');
sex =4; % chose sex 3 4 5 for female male total
x_max=102; %new maximum age group 100+, set to prevent 0 and NaN in M
x_min= 0; %new mimimum age group, set to 0 skip and use all grous 0-x_max
t_min=1956; %first year of data included
t_max=2015; %last year of data included
UseAccumulatedEndGroup=0; %1 yes 0 no to calc new accumulated m(x,t) for group x_max+ 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matrix m(x,t) dimensions
x=x_max-x_min+1; %number of age groups used
t=t_max-t_min+1; % number of years used

%autodetected x_full for maximal age in data, i.e 110+
%autodetected t for age interval of data

DeathRateData=[DeathRateData.colheaders;num2cell(DeathRateData.data)];
Death_Matrix=cell2mat(DeathRateData(2:end,:));

t_start=cell2mat(DeathRateData(2,1)); %data first recored year, 1958
t_end=cell2mat(DeathRateData(end,1)); %data last recorded year, 2015
t_full=t_end-t_start+1; %data total number of years recorded, 60

%x_full=110+1; maxmium number of age groups ins data
x_full=cell2mat(DeathRateData(end,2))+1; %maximum age group 110+, i.e. x_full-1+

%Central Death Matrix m(x,t)
M_full=zeros(x_full,t_full);

for i = 1:t_full
   M_full(:,i)=Death_Matrix(1+(i-1)*x_full:x_full*i,sex);  
end

%for weighted fit we replace -1 and 0 with an arbitrary number
if(replace42==1)
    B=M_full<=0;
    M_full(B)=42;
end

M=zeros(x,t);
M(1:x,1:t)=M_full(x_min+1:x_max+1,(t_min-t_start+1):t_full-(t_end-t_max));

%central death rate M(x,t)= D(x,t)/E(x,t)
%deaths/exposure to risk
%M(100+,t)=Sum D(100+,t)/ Sum E(100+,t) = CDM(101,t)

E=[E.colheaders;num2cell(E.data)];
Exposures=cell2mat(E(2:end,:));
D=[D.colheaders;num2cell(D.data)];
Deaths=cell2mat(D(2:end,:));

%D_xt E_xt matrixes with group 100+ containing all info of groups 100-110+
D_full=zeros(x_full,t_full);
E_full=zeros(x_full,t_full);

for i = 1:t_full
   D_full(:,i)=Deaths(1+(i-1)*x_full:x_full*i,sex);  
   E_full(:,i)=Exposures(1+(i-1)*x_full:x_full*i,sex);  
end
D=zeros(x,t);
D(1:x,1:t)=D_full(x_min+1:x_max+1,(t_min-t_start+1):t_full-(t_end-t_max));
E=zeros(x,t);
E(1:x,1:t)=E_full(x_min+1:x_max+1,(t_min-t_start+1):t_full-(t_end-t_max));

%zeros and missing values for ages above 100 make it hard to work with 
%If CalcAccumulatedEndGroup=1 recalculate from the population data 0-110+ the central death rate for a new group
%100+ and use this for further calculations 
%else just cut of the data at age 100
if UseAccumulatedEndGroup==1
    for i=1:t %years from t_min until t_start, i.e. 1956-2015 for full data
        D(x,i)=sum(D_full(x_max+1:end,i+t_min-t_start));
        E(x,i)=sum(E_full(x_max+1:end,i+t_min-t_start));
        M(x,i)=D(x,i)/E(x,i);

    end
end
lnM=log(M);

clearvars i j DeathRateData Exposures Deaths Death_Matrix B
%clearvars t_start t_end t_full x_full M_full 
