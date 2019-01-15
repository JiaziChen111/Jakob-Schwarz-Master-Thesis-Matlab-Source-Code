xx=linspace(1,x,x);
tt=linspace(1,t,t);

figure;

subplot(3,2,1);
plot(xx,a_wLC,'r')
hold on;
plot(xx,a_gdp_normalen,'c');
hold on;
plot(xx,a_gdp_newton,'b');
hold on;
plot(xx,a_gdp_simplex,'g');
hold on;
plot(xx,a_gdp_newton2,'k');
hold on;
%plot(xx,a_gdp_simplex2,'y');

title('a_x');
legend('wLC','normalen','newton','simplex','newton2');
%legend('wLC','newton','simplex','newton2','simplex2');

subplot(3,2,2);
plot(xx,b_wLC,'r')
hold on;
plot(xx,b_gdp_normalen,'c');
hold on;
plot(xx,b_gdp_newton,'b');
hold on;
plot(xx,b_gdp_simplex,'g');
hold on;
plot(xx,b_gdp_newton2,'k');
hold on;
%plot(xx,b_gdp_simplex2,'y');

title('b_x');
legend('wLC','normalen','newton','simplex','newton2');
%legend('wLC','newton','simplex','newton2','simplex2');

subplot(3,2,3);
plot(tt,k_wLC,'r')
hold on;
plot(tt,k_gdp_normalen,'c');
hold on;
plot(tt,k_gdp_newton,'b');
hold on;
plot(tt,k_gdp_simplex,'g');
hold on;
plot(tt,k_gdp_newton2,'k');
hold on;
%plot(tt,k_gdp_simplex2,'y');
title('k_t');

legend('wLC','normalen','newton','simplex','newton2');
%legend('wLC','newton','simplex','newton2','simplex2');


subplot(3,2,4);
plot(xx,y_gdp_normalen,'c');
hold on;
plot(xx,y_gdp_newton,'b');
hold on;
plot(xx,y_gdp_simplex,'g');
hold on;
plot(xx,y_gdp_newton2,'k');
hold on;
%plot(xx,y_gdp_simplex2,'y');

title('y_x');
legend('normalen','newton','simplex','newton2');
%legend('newton','simplex','newton2','simplex2');


subplot(3,2,5);
plot(tt,g_t);
title('g_t');

clearvars xx tt