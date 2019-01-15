function [ e_xt, e_xt_UB, e_xt_LB ] = lifeexp_m(X, T, M, x_min, x_max, Msim,alpha)
%LIFEEXP Calculates the remaining full years life expectancy
%   X is a vector containing current ages at times T 
%   X = [1:x] entspricht [x_min:x_max] used in q
%   T is a vector of times [1:t] entspricht [t_min:t_max] of q
%   e_xt is the matrix of all combinations of X and T
%   M is the central death rate matrix q(x,t)
%   OPTINAL Msim contains a cell array of M used to get CI
%   Formula from Poison paper, max k=1:x_max-x (such that k+x=x_max)
e_xt=zeros(length(X),length(T));
for i=1:length(X)
   %actual we start with
   x=x_min-1+X(i);
   for j=1:length(T)
      %prob to survive 1 year aged X(i) year T(i)
      k_p_xt=1;
      e_xt(i,j)=(1-exp(-M(X(i),T(j))))/M(X(i),T(j));
      for k=1:x_max-x
        %prob to survive k-1 years*prob to survive 1 year 
        k_p_xt=k_p_xt*exp(-M(X(i)+k-1,T(j)+k-1));
        e_xt(i,j)=e_xt(i,j)+k_p_xt * (1-exp(-M(X(i)+k,T(j)+k)))/M(X(i)+k,T(j)+k);
        
      end
   end
    
end
%test if optinal Msim argument is given, if yes calculate from it CI
if nargin==7
   e_xt_UB=zeros(length(X),length(T));
   e_xt_LB=zeros(length(X),length(T));
   H=zeros(length(Msim),1);
   for i=1:length(X)
   %actual we start with
   x=x_min-1+X(i);
       for j=1:length(T)
           %for each patch calculate e_xt and safe in H
           for n=1:length(Msim)
           
              %prob to survive 1 year aged X(i) year T(i)
              k_p_xt=1;
              H(n)=(1-exp(-Msim{n}(X(i),T(j))))/Msim{n}(X(i),T(j));
              for k=1:x_max-x
                %prob to survive k-1 years*prob to survive 1 year 
                k_p_xt=k_p_xt*exp(-Msim{n}(X(i)+k-1,T(j)+k-1));
                H(n)=H(n)+k_p_xt * (1-exp(-Msim{n}(X(i)+k,T(j)+k)))/Msim{n}(X(i)+k,T(j)+k);
              end   
           end
       e_xt_LB(i,j)=prctile(H',100*(alpha/2),2);  
       e_xt_UB(i,j)=prctile(H',100*(1-alpha/2),2);
       end
   end  
end


end

