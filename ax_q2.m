function [ e_xt, e_xt_UB, e_xt_LB ] = ax_q2(X, T, q_M, x_min, x_max, q_future_mean,q_future_paths,alpha,irate)
%LIFEEXP Calculates the remaining full years life expectancy
%   X is a vector containing current ages at times T 
%   X = [1:x] entspricht [x_min:x_max] used in q
%   T is a vector of times [1:t] entspricht [t_min:t_max] of q
%   e_xt is the matrix of all combinations of X and T
%   q is the death prob matrix q(x,t)
%   Sum (k=0:x_max-x) k* [k_p_x(t)]*[q_x+k(t+k)]   + [(z-x+1)* (x_max+1)_p_x(t)*q=1],
%   OPTIONAL if A_future given with paths then a CI is calculated
e_xt=ones(length(X),length(T));
q=[q_M,q_future_mean];
n=length(q_future_paths);

for i=1:length(X)
   %actual we start with
   x=x_min-1+X(i);
   for j=1:length(T)
      %prob to survive 1 year aged X(i) year T(i)
      k_p_xt=1-q(X(i),T(j));
      vk=1/(1+irate);
      for k=1:x_max-x
        e_xt(i,j)=e_xt(i,j)+vk*k_p_xt;
        %prob to survive k years *[prob to survive to next year]
        k_p_xt=k_p_xt*(1-q(X(i)+k,T(j)+k));
        vk=vk*vk;
      end
      %for example 0-100 then prob to reach 101 is positive and known from 
      %data and therefor can need to be included
      %if age is 101 (or above) prob to reach 102 is zero
      if x<=x_max
         e_xt(i,j)=e_xt(i,j)+vk*k_p_xt;
      end
   end   
end

%test if optinal Msim argument is given, if yes calculate from it CI
   
   e_xt_UB=zeros(length(X),length(T));
   e_xt_LB=zeros(length(X),length(T));
   
   for i=1:length(X)
   %actual we start with
   x=x_min-1+X(i);
       for j=1:length(T)
           H=ones(n,1);
           %for each patch calculate e_xt and safe in H
           for r=1:n 
               vk=1/(1+irate);
               qsim=[q_M,q_future_paths{r}]; 
               %prob to survive 1 year aged X(i) year T(i)
              k_p_xt=1-qsim(X(i),T(j));
              for k=1:x_max-x
                H(r)=H(r)+vk*k_p_xt;
                %prob to survive k years *[prob to survive to next year]
                k_p_xt=k_p_xt*(1-qsim(X(i)+k,T(j)+k));
                vk=vk*vk;
              end
              %for example 0-100 then prob to reach 101 is positive, dying prob set to 1
              %if age is 101 (or above) prob to reach 102 is zero
              if x<=x_max
                 H(r)=H(r)+vk*k_p_xt;
              end   
           end
       e_xt_LB(i,j)=prctile(H',100*(alpha/2),2);  
       e_xt_UB(i,j)=prctile(H',100*(1-alpha/2),2);
       end
   end  
end

