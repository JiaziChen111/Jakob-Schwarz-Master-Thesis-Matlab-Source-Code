function [ e_xt, e_xt_UB, e_xt_LB , qsim] = lifeexp_q(X, T, q, x_min, x_max, A_future,alpha)
%LIFEEXP Calculates the remaining full years life expectancy
%   X is a vector containing current ages at times T 
%   X = [1:x] entspricht [x_min:x_max] used in q
%   T is a vector of times [1:t] entspricht [t_min:t_max] of q
%   e_xt is the matrix of all combinations of X and T
%   q is the death prob matrix q(x,t)
%   Sum (k=0:x_max-x) k* [k_p_x(t)]*[q_x+k(t+k)]   + [(z-x+1)* (x_max+1)_p_x(t)*q=1],
%   OPTIONAL if A_future given with paths then a CI is calculated
e_xt=zeros(length(X),length(T));
for i=1:length(X)
   %actual we start with
   x=x_min-1+X(i);
   for j=1:length(T)
      %prob to survive 1 year aged X(i) year T(i)
      k_p_xt=1-q(X(i),T(j));
      for k=1:x_max-x
        e_xt(i,j)=e_xt(i,j)+k*k_p_xt * q(X(i)+k,T(j)+k);
        %prob to survive k years *[prob to survive to next year]
        k_p_xt=k_p_xt*(1-q(X(i)+k,T(j)+k));
      end
      %for example 0-100 then prob to reach 101 is positive, dying prob set to 1
      %if age is 101 (or above) prob to reach 102 is zero
      if x<=x_max
         e_xt(i,j)=e_xt(i,j)+(x_max+1-x)*k_p_xt;
      end
   end   
end

%test if optinal Msim argument is given, if yes calculate from it CI
if nargin==7
    t_future=length(A_future{1}); %number of future years
    h=size(q);
    %cut out only the fit matrix
    q=q(:,1:h(2)-t_future);
    %%%calculate Qsim cell array with paths from A_future
    
    n=length(A_future); %number of paths
    qsim=cell(n,1); %each cell is a q path
    x=size(q); 
    t_C=x(2); %ages for fit
    x=x(1); %ages
    
    
    q_Forecast=zeros(x,t_future);
    %outer path loop
    for k=1:n
        for i=1:x
            %forecasted year loop
            for j=1:t_future
                h=A_future{k}(j,1)+A_future{k}(j,2)*((x_min+i-1)+(j-1+t_C));
                q_Forecast(i,j)=exp(h)/(1+exp(h));
            end
        end
        qsim{k}=[q,q_Forecast]; 
    end
   %%%
   
   e_xt_UB=zeros(length(X),length(T));
   e_xt_LB=zeros(length(X),length(T));
   
   for i=1:length(X)
   %actual we start with
   x=x_min-1+X(i);
       for j=1:length(T)
           H=zeros(n,1);
           %for each patch calculate e_xt and safe in H
           for r=1:n  
               %prob to survive 1 year aged X(i) year T(i)
              k_p_xt=1-qsim{r}(X(i),T(j));
              for k=1:x_max-x
                H(r)=H(r)+k*k_p_xt * qsim{r}(X(i)+k,T(j)+k);
                %prob to survive k years *[prob to survive to next year]
                k_p_xt=k_p_xt*(1-qsim{r}(X(i)+k,T(j)+k));
              end
              %for example 0-100 then prob to reach 101 is positive, dying prob set to 1
              %if age is 101 (or above) prob to reach 102 is zero
              if x<=x_max
                 H(r)=H(r)+(x_max+1-x)*k_p_xt;
              end   
           end
       e_xt_LB(i,j)=prctile(H',100*(alpha/2),2);  
       e_xt_UB(i,j)=prctile(H',100*(1-alpha/2),2);
       end
   end  
end

end

