function [a_x,b_x,k_t,y_x,fval]=PoissonFit_explanation(D, E,g_t, a_start,b_start,k_start,y_start,norming)
%PoissonFit calculates a Maximum likelihood estimation by a uni dimensional
%elemntary Newton method given by "A Poisson log-bilinear regression approach to
%the construction of projected lifetables by Natacha Brouhns, Michel Denuit,Jeroen K. Vermunt"
%For their model D_xt~Poi(E_xt* m(x,t)) with m(x,t)=e^(a_x+b_x*k_t)
%max L(a,b,k)=Sum(x,t) D_xt (a_x+b_x*k_t)- E_xt exp(a_x+b_x*k_) + constant

 %g_t is a matrix= [g_1(1:n)|...|g_t(1:n)]
    %y_x is a matrix= [y_1:x(1)|...|y_1:x(n)]
    %norming=1 cov(k_t-k_t-1,g_t-g_t-1)=0  else cov(k_t,g_t)=0

%starting values
n=size(D);
x=n(1);
t=n(2);
n=size(y_start);
n=n(2);

   a_x=a_start;
   b_x=b_start;
   k_t=k_start;
   y_x=y_start;
   Hg=y_x*g_t;

%log-likelihood funtion to maximise used as break condition
     function L = L(a_x,b_x,k_t)
        L=0;
        for ix=1:x
            for j=1:t
                L=L+D(ix,j)*(a_x(ix)+b_x(ix)*k_t(j)+Hg(ix,j))-E(ix,j)*exp(a_x(ix)+b_x(ix)*k_t(j)+Hg(ix,j) );
            end
        end
     end
L_new=L(a_start,b_start,k_start);
L_old=L_new-1000;

%help variables
h=zeros(x,1);
hh=zeros(t,1);
hhh=zeros(x,n);
%it is possible that we never start the while loop due to L_old choice
while (abs(L_new-L_old)>10e-6) %look if change of L is small
    %estimate D_hat
    D_hat=zeros(x,t);
        for ii=1:x
            for jj=1:t
                D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj)+Hg(ii,jj) );
            end
        end
    Diff=D-D_hat;
    %update a_x  
    for iii=1:x
        h(iii)=sum (Diff(iii,:))/(-sum(D_hat(iii,:)));
    end
    a_x=a_x-h;
    %update estimate D_hat
    D_hat=zeros(x,t);
        for ii=1:x
            for jj=1:t
                 D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj)+Hg(ii,jj) );
            end
        end
    Diff=D-D_hat;
    %update k
    for jjj=1:t
        hh(jjj)=sum(Diff(:,jjj).*b_x)/(-sum(D_hat(:,jjj).*(b_x.^2)));
    end
    k_t=k_t-hh; 
    %update estimate D_hat
    D_hat=zeros(x,t);
        for ii=1:x
            for jj=1:t
                 D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj)+Hg(ii,jj) );
            end
        end
    Diff=D-D_hat;
    %update b
    for iii=1:x
        h(iii)=sum(Diff(iii,:).*k_t')/(-sum(D_hat(iii,:).*(k_t.^2)'));
    end
    b_x=b_x-h;    
    
    %update estimate D_hat
    D_hat=zeros(x,t);
        for ii=1:x
            for jj=1:t
                 D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj)+Hg(ii,jj) );
            end
        end
    Diff=D-D_hat;
    %update y
    
    for k=1:n
        for iii=1:x
            hhh(iii,k)=sum(Diff(iii,:).*g_t(k,:))/(-sum(D_hat(iii,:).*(g_t(k,:).^2)));
        end
    end
    y_x=y_x-hhh;    
    
    %update Hg
    Hg=y_x*g_t;
       
    %update break condition
    L_old=L_new;
    L_new=L(a_x,b_x,k_t); 
end
fval=L(a_x,b_x,k_t);

%norming
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

