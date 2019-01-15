function [ a_x,b_x,k_t ] = LeeCarterFit(D, E, a, b, k_start)
%SVDFit adjust the SVD Fits' k_t to reproduce the observed number of deaths
%   The following is losely cited from "A Poisson log-bilinear regression approach to
%   the construction of projected lifetables by Natacha Brouhns, Michel
%   Denuit and Jeroen K. Vermunt":

%   Solve Sum(x) D_xt= Sum(x) E_xt* e^(a_x+b_x*k_t) for k_t
%   resulting death rates produce the total number of deaths actually observed 
%   in the data for the year t in question
%   There are several advantages to make this second-stage estimate of the
%   parameters k_t
%   In particular, it avoids sizable discrepancies between predicted and actual deaths 
%   (occurring because the first step is based on logarithms of death
%   rates). Other advantages are discussed by Lee (2000).

%start value for linesearch
k_t=k_start;
x=length(a);     
     %error funtion to minimise
     function err = Errfun(k,t) %error function for min search of k_t,
        %t=1-60: 1956-2015
        err=0; % weighted sum of all exposure year t
        for j=1:x
            err=err+E(j,t)*exp((a(j)+b(j)*k));
        end
        err=abs(sum (D(:,t))-err); %sum of all deaths year t-SumE);
     end
    
 for t=1:length(k_start)
     FUN = (@(k) Errfun(k,t));     
     %line search with start value k by Nelder-Meat
     k_t(t) = fminsearch(FUN,k_start(t));
     
 end

%normalize such that sum b =1
c=sum(b);
b_x=b/c;
k_t=k_t*c;

%sum(k_t)=0
c=-sum(k_t)/t;
a_x=a-b*c;
k_t=k_t+c;

 
end


