function  [a_x,b_x,k_t,fval,exitflag,steps]=Fit_wLC(lnM,D)
%FIT_WLC fit the Lee Carter model by weighted least sqarres
%   a solution that minimises the quadratic error is first calculated
%   D is optional for weighted least squarres. if not given set to ones
%   afterwards the solution is transfomred into a nomalized solution
%   introducing weighted least sqarres erases zero cell problems
%   min sum(x,t) w_xt (ln m(x,t)-a_x-b_x k_t      (1)
%   w_xt=D_xt number of deaths

%   to minimise (1) calculate first derivates w.r.t a_x, b_x, k_t
%   yields normal euqations for fixed x bzw t: a_x=...,b_x=.. k_t=...
%   solve while parameter change is big
        % set a_x=...for all x
        % set b_x=...for all x
        % set k_t=...for all t
    %end
%alternatively solve (1) with quasi newton or simplex
%Wilmoth 1993 found quasi newton to be fastet and simplex very slow
%solution error is almost equal for all three with normal eqations best
%then simplex and worst quasi newton.
n=size(lnM);
x=n(1);
t=n(2);

%if no weights given set to 1 except for where 42 (M<=0)
if nargin==1
    D=ones(n);
    B=lnM==log(42);
    D(B)=0;
end

a_x=zeros(x,1);
b_x=ones(x,1);
k_t=ones(t,1);

a_temp=ones(x,1);
b_temp=zeros(x,1);
k_temp=ones(t,1);
exitflag=1;
steps=1;
%minimum amount of steps can be set
while(norm(a_x-a_temp)>10e-6 || norm(b_x-b_temp)>10e-6 || norm(k_t-k_temp)>10e-6 || steps<3)
    a_temp=a_x;
    b_temp=b_x;
    k_temp=k_t;
    
    for ii=1:x
        a_x(ii)=sum(D(ii,:).*(lnM(ii,:)-b_x(ii)*k_t'))/sum(D(ii,:));
    end
    
    for ii=1:x
        b_x(ii)=sum(D(ii,:).*k_t'.*(lnM(ii,:)-a_x(ii)))/sum(D(ii,:).*k_t'.^2);
    end
    
    for jj=1:t
         k_t(jj)=sum(D(:,jj).*b_x.*(lnM(:,jj)-a_x))/sum(D(:,jj).*b_x.^2);
    end
    steps=steps+1;   
     h=sum(a_x)+sum(b_x)+sum(k_t);
    if  isnan(h) || isinf(h)
    a_x=a_temp;
    b_x=b_temp;
    k_t=k_temp;
    exitflag=-1;
    end
end


%normalize such that sum b =1
c=sum(b_x);
b_x=b_x/c;
k_t=k_t*c;

%sum(k_t)=0
c=-sum(k_t)/t;
a_x=a_x-b_x*c;
k_t=k_t+c;

FUN = @(z) sum(sum(D.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)').^2   ))   ; 
x0 = [a_x;b_x;k_t];
fval=FUN(x0);
end

