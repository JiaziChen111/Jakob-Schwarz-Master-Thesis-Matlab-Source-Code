function [a_x,b_x,k_t,y_x,fval,exitflag_normalen,steps]=FitLC_GDP_normalen(lnM,g_t, a_start, b_start, k_start,y_start,D)
%FITLC_GDO fit the extended GDP Lee Carter model: ln m(x,t)=a_x+b_x*k_t+y_x*g_t+e
    %D optinal weights    
    %see Fit_wLC.m for more details.


n=size(lnM);
x=n(1);
t=n(2);

a_x=a_start;
b_x=b_start;
k_t=k_start;
y_x=y_start;


    function [a,b,k,y]= normalize(a,b,k,y)
        e=cov(k,g_t);
        e=e(1,2)/e(2,2);
        d=sum(b);
        k=d*(k-e*g_t);
        y=y+e*b;
        b=b/d;

        c=-sum(k)/t;
        a=a-b*c;
        k=k+c;
    end

%[a_x,b_x,k_t,y_x]= normalize(a_x,b_x,k_t,y_x);

 a_temp=a_x+ones(x,1);
 b_temp=b_x;
 k_temp=k_t;
 y_temp=y_x;
 steps=1;
 exitflag_normalen=1;
while(norm(a_x-a_temp)>10e-10 || norm(b_x-b_temp)>10e-15|| norm(k_t-k_temp)>10e-10 ||  norm(y_x-y_temp)>10e-10)
    a_temp=a_x;
    b_temp=b_x;
    k_temp=k_t;
    y_temp=y_x;
    
    
    for ii=1:x
        a_x(ii)=sum(D(ii,:).*(lnM(ii,:)-b_x(ii)*k_t'-y_x(ii)*g_t'))/sum(D(ii,:));
    end
    
    for ii=1:x
        b_x(ii)=sum(D(ii,:).*k_t'.*(lnM(ii,:)-a_x(ii))-y_x(ii)*g_t')/sum(D(ii,:).*k_t'.^2);
    end
    
    for jj=1:t
        k_t(jj)=sum(D(:,jj).*b_x.*(lnM(:,jj)-a_x-y_x*g_t(jj)))/sum(D(:,jj).*b_x.^2);
    end
     [a_x,b_x,k_t,y_x]= normalize(a_x,b_x,k_t,y_x);
    for ii=1:x
        y_x(ii)=sum(D(ii,:).*g_t'.*(lnM(ii,:)-a_x(ii))-b_x(ii)*k_t')/sum(D(ii,:).*g_t'.^2);
    end
   
  
    steps=steps+1; 
    
    h=sum(a_x)+sum(b_x)+sum(k_t)+sum(y_x);
    if  isnan(h) || isinf(h)
    a_x=a_temp;
    b_x=b_temp;
    k_t=k_temp;
    y_x=y_temp; 
    exitflag_normalen=-1;
    end
    
end

[a_x,b_x,k_t,y_x]= normalize(a_x,b_x,k_t,y_x);
FUN = @(z) sum(sum(D.*(lnM-repmat(z(1:x),1,t)-z(x+1:2*x)*z(2*x+1:2*x+t)'-z(2*x+t+1:end)*g_t').^2   ))   ; 
fval=FUN([a_x;b_x;k_t;y_x]);

end

