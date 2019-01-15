function [a_x,b_x,k_t,y_x] = explanation_renorm(a_x,b_x,k_t,g_t ,norming)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%norming
x=length(a_x);
t=length(k_t);
n=size(g_t);
n=n(1);
y_x=zeros(x,n);
if (norming ==1) %differene independece neu
     dg=diff(g_t');
     CovMatrixdgti= cov(dg);
     covVectorktdgti=zeros(1,n);
     for i=1:n
         temp=cov(diff(k_t),diff(g_t(i,:)));
         covVectorktdgti(i)=temp(2,1);      
     end
         
     e=CovMatrixdgti\covVectorktdgti;      
     
else %normal alt
     CovMatrixgti= cov(g_t');
     covVectorktgti=zeros(n,1);
     for i=1:n
         temp=cov(k_t,g_t(i,:));
         covVectorktgti(i)=temp(2,1);      
     end
     e=CovMatrixgti\covVectorktgti;     
    
end
        
        k_t=k_t-g_t'*e;
        y_x=y_x+b_x*e';     
        h=sum(b_x);
        b_x=b_x/h;
        k_t=k_t*h;
        c=-sum(k_t)/t;
        a_x=a_x-b_x*c;
        k_t=k_t+c;        
    
end
