function [a_x,b_x,k_t,y_x ] = LeeCarterFit_Explanation(D, E, a_x, b_x, k_start, y_x, g, norming)
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
x=length(a_x); 
n=size(g);
n=n(1);
     %error funtion to minimise
     function err = Errfun(k,t) %error function for min search of k_t,
        %t=1-60: 1956-2015
        err=0; % weighted sum of all exposure year t
        for j=1:x
            err=err+E(j,t)*exp((a_x(j)+b_x(j)*k)+y_x(j,:)*g(:,t));
        end
        err=abs(sum (D(:,t))-err); %sum of all deaths year t-SumE);
     end
    
 for t=1:length(k_start)
     FUN = (@(k) Errfun(k,t));     
     %line search with start value k by Nelder-Meat
     k_t(t) = fminsearch(FUN,k_start(t));
     
 end

%norming
if (norming ==1) %my new differene independece
     dg=diff(g');
     CovMatrixdgti= cov(dg);
     covVectorktdgti=zeros(1,n);
     for i=1:n
         temp=cov(diff(k_t),diff(g(i,:)));
         covVectorktdgti(i)=temp(2,1);      
     end
         
     e=CovMatrixdgti\covVectorktdgti;      
     
else %normal standard not mine
     CovMatrixgti= cov(g');
     covVectorktgti=zeros(n,1);
     for i=1:n
         temp=cov(k_t,g(i,:));
         covVectorktgti(i)=temp(2,1);      
     end
     e=CovMatrixgti\covVectorktgti;     
    
end     
        k_t=k_t-g'*e;
        y_x=y_x+b_x*e';     
        h=sum(b_x);
        b_x=b_x/h;
        k_t=k_t*h;
        c=-sum(k_t)/t;
        a_x=a_x-b_x*c;
        k_t=k_t+c;    
 
end


