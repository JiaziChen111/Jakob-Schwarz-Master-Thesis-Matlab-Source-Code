%changes to adapt for other data see Master for setting:
%{
x_max_C=109; %new maximum age considered, use 109 to exlude last row of 1
x_min_C= 40; %new mimimum age considered, Cairns starts with 60, 40 seesm to work well already
t_min_C=1956; %first year of data included
t_max_C=2015; %last year of data included
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sex= 3, 4, 5 for female, male, total, choose in BuildCentralDeathMatrix
switch sex
    case 3
        LifetableData= importdata('fltper_1x1.txt');
    case 4
        LifetableData= importdata('mltper_1x1.txt');
    case 5
        LifetableData= importdata('bltper_1x1.txt');
    otherwise
        disp('invalid gender choice for variable sex')
end

LifetableData=[LifetableData.colheaders;num2cell(LifetableData.data)];
Lifetable=cell2mat(LifetableData(2:end,:));

%matrix q(x,t) dimensions
x_C=x_max_C-x_min_C+1; %number of age groups used
t_C=t_max_C-t_min_C+1; % number of years used

t_start_C=Lifetable(1,1); %data first recored year, 1958
t_end_C=Lifetable(end,1); %data last recorded year, 2015
t_full_C=t_end-t_start_C+1; %data total number of years recorded, 60

%x_full=110+1; maxmium number of age groups ins data
x_full_C=Lifetable(end,2)+1; %maximum age group 110+, i.e. x_full-1+

%1 year Death Probability Matrix q(x,t)
q_full=zeros(x_full_C,t_full_C);
for i = 1:t_full_C
   q_full(:,i)=Lifetable(1+(i-1)*x_full_C:x_full_C*i,sex);  
end

%Cairns only used death prob of age 60+
q=zeros(x_C,t_C);
q(1:x_C,1:t_C)=q_full(x_min_C+1:x_max_C+1,(t_min_C-t_start_C+1):t_full_C-(t_end_C-t_max_C));

clearvars LifetableData Lifetable i 
%clearvars t_end_C t_full_C x_full_C t_start_C q_full