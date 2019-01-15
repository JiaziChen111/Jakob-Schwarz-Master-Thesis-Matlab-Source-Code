function [a_x,b_x,k_t,fval]=PoissonFit(D, E, a_start,b_start,k_start)
%PoissonFit calculates a Maximum likelihood estimation by a uni dimensional
%elemntary Newton method given by "A Poisson log-bilinear regression approach to
%the construction of projected lifetables by Natacha Brouhns, Michel Denuit,Jeroen K. Vermunt"
%For their model D_xt~Poi(E_xt* m(x,t)) with m(x,t)=e^(a_x+b_x*k_t)
%max L(a,b,k)=Sum(x,t) D_xt (a_x+b_x*k_t)- E_xt exp(a_x+b_x*k_) + constant

%starting values
n=size(D);
x=n(1);
t=n(2);

if nargin==5 
    a_x=a_start;
    b_x=b_start;
    k_t=k_start;
else
    a_x=zeros(x,1);
    b_x=ones(x,1);
    k_t=zeros(t,1);
end

%log-likelihood funtion to maximise used as break condition
     function L = L(a_x,b_x,k_t)
        L=0;
        for i=1:x
            for j=1:t
                L=L+D(i,j)*(a_x(i)+b_x(i)*k_t(j))-E(i,j)*exp(a_x(i)+b_x(i)*k_t(j));
            end
        end
     end
L_new=L(a_x,b_x,k_t);
L_old=L_new-1000;

%help variables
h=zeros(x,1);
hh=zeros(t,1);
%it is possible that we never start the while loop due to L_old choice
while (abs(L_new-L_old)>10e-10) %look if change of L is small
    %estimate D_hat
    D_hat=zeros(x,t);
        for ii=1:x
            for jj=1:t
                D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj));
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
                D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj));
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
                D_hat(ii,jj)=E(ii,jj)*exp(a_x(ii)+b_x(ii)*k_t(jj));
            end
        end
    Diff=D-D_hat;
    %update b
    for iii=1:x
        h(iii)=sum(Diff(iii,:).*k_t')/(-sum(D_hat(iii,:).*(k_t.^2)'));
    end
    b_x=b_x-h;    
    %update break condition
    L_old=L_new;
    L_new=L(a_x,b_x,k_t); 
end
fval=L(a_x,b_x,k_t);
%following part is a improvisation of me and not from the original paper

%ensuring model indetification sum(x) b_x = 1 
%note for any c is a,b*c,k/c also a solution for ln m_xt=a+b*k
c=sum(b_x);
b_x=b_x/c;
k_t=k_t*c;
%ensuring model indetification sum k=0
%note a -bc, b , k+c is also a solution
c=sum(k_t)/t;
k_t=k_t-c;
a_x=a_x+b_x*c;

end

