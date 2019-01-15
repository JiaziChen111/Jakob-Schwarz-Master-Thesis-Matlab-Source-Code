function [a_x,b_x,k_t,y_x,fval,exitflag,steps] = FitLC_Explanation_2step(lnM,g,norming,D,E)
%FitLC_Explanation_2step Uses 2-step to find a solution to the LC
%explanation model either input W or D E
%if D is weights used for wLC in second step
% D E for ML fit, if those supplied then ML fit else wLC fit
n=size(lnM);
x=n(1);
t=n(2);
n=size(g);
n=n(1);

tilde_a_x=zeros(x,1);
%b_x=ones(x,1);
%k_t=ones(t,1);
y_x=zeros(x,n);
index42=lnM==log(42);
index42=logical((index42-1)*(-1));
%step 1
for i=1:x
    H=index42(i,:);
    Y=lnM(i,:)';
    X=[ones(t,1),g'];  
    %H is used to exlude zero cells
    b=X(H,:)\Y(H); 
    tilde_a_x(i)=b(1);
    y_x(i,:)=b(2:end);
end
%step 2
lnM2=lnM-repmat(tilde_a_x,1,t);
for i=1:n
    lnM2=lnM2-y_x(:,i)*g(i,:);    
end
%dont forget zero cells
index42=lnM==log(42);
lnM2(index42)=log(42);

%%% solve Lee Carter system
if nargin==4 %wLC with weights D
    [a_x,b_x,k_t,~,exitflag,steps]=Fit_wLC(lnM2,D);
    a_x=a_x+tilde_a_x; %step 3
        %%%calculate fitting error
        FUN = @(z) sum(sum(D.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)'-reshape(z(2*x+t+1:end),x,n)*g).^2   ));
        y0=reshape(y_x,x*n,1);
        x0 = [a_x;b_x;k_t;y0];
        fval=FUN(x0);

        %ML does not take into acount lnM2
%elseif nargin==5 %ML with D and E 
%    [a_x,b_x,k_t,~]=PoissonFit(D, E);
%    a_x=a_x+tilde_a_x; %step 3
%    L = @(z) sum(sum( E.*exp(repmat(z(1:x),1,t)).*exp(z(x+1:2*x)*z(2*x+1:2*x+t)').*exp(reshape(z(2*x+t+1:end),x,n)*g_t) - D.*(repmat(z(1:x),1,t)+z(x+1:2*x)*z(2*x+1:2*x+t)'+reshape(z(2*x+t+1:end),x,n)*g_t)  ))   ; 
%    y0=reshape(y_start,x*n,1);   
%    x0 = [a_start;b_start;k_start;y0];
%    fval=L(x0);
%    fval=-fval;

else %wLC with 1-0 weights  
        D=ones(x,t);
        B=lnM==log(42);
        D(B)=0;
    [a_x,b_x,k_t,~,exitflag,steps]=Fit_wLC(lnM2,D);
    a_x=a_x+tilde_a_x; %step 3
        %%%calculate
        FUN = @(z) sum(sum(D.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)'-reshape(z(2*x+t+1:end),x,n)*g).^2   ));
        y0=reshape(y_x,x*n,1);
        x0 = [a_x;b_x;k_t;y0];
        fval=FUN(x0);
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

