figure;
plot(A_1,A_2);
xlabel('A_1');
ylabel('A_2');
title('A(t))');

figure;
dA=[diff(A_1),diff(A_2)];
plot(dA(:,1),dA(:,2),'*');
title('A(t)-A(t-1)');

figure;
plot(linspace(1,x,x),a_LC,'r');
hold on;
plot(linspace(1,x,x),a_SVD,'k');
hold on;
plot(linspace(1,x,x),a_ML,'b');
hold on;
plot(linspace(1,x,x),a_wLC,'c');
hold on;
plot(linspace(1,x,x),a_wLC2,'m');
title('a_x');
legend('LC','SVD','ML','wLC','wLC2');

figure;
plot(linspace(1,x,x),b_LC,'r');
hold on;
plot(linspace(1,x,x),b_SVD,'k');
hold on;
plot(linspace(1,x,x),b_ML,'b');
hold on;
plot(linspace(1,x,x),b_wLC,'c');
hold on;
plot(linspace(1,x,x),b_wLC2,'m');
title('b_x');
legend('LC','SVD','ML','wLC','wLC2');

figure('name','ACD and PACF for k_t');
subplot(6,2,[1 2]);
plot(linspace(1,t,t),k_LC,'r');
hold on;
plot(linspace(1,t,t),k_SVD,'k');
hold on;
plot(linspace(1,t,t),k_ML,'b');
hold on;
plot(linspace(1,t,t),k_wLC,'c');
hold on;
plot(linspace(1,t,t),k_wLC2,'m');
title('k_t');
legend('LC','SVD','ML','wLC','wLC2');

subplot(6,2,3) 
autocorr(k_LC)
title('Sample Autocorrelation Function for k_{LC}')
subplot(6,2,4)
parcorr(k_LC)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for k_{LC}')

subplot(6,2,5) 
autocorr(k_SVD)
title('Sample Autocorrelation Function for k_{SVD}')
subplot(6,2,6)
parcorr(k_SVD)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for k_{SVD}')

subplot(6,2,7) 
autocorr(k_ML)
title('Sample Autocorrelation Function for k_{ML}')
subplot(6,2,8)
parcorr(k_ML)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for k_{ML}')

subplot(6,2,9) 
autocorr(k_wLC)
title('Sample Autocorrelation Function for k_{wLC}')
subplot(6,2,10)
parcorr(k_wLC)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for k_{wLC}')

subplot(6,2,11) 
autocorr(k_wLC2)
title('Sample Autocorrelation Function for k_{wLC2}')
subplot(6,2,12)
parcorr(k_wLC2)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for k_{wLC2}')


figure('name','ACD and PACF for dk_t');
dk_LC=diff(k_LC);
dk_SVD=diff(k_SVD);
dk_ML=diff(k_ML);
dk_wLC=diff(k_wLC);
dk_wLC2=diff(k_wLC2);
subplot(6,2,[1 2]);
plot(linspace(1,t-1,t-1),dk_LC,'r');
hold on;
plot(linspace(1,t-1,t-1),dk_SVD,'k');
hold on;
plot(linspace(1,t-1,t-1),dk_ML,'b');
hold on;
plot(linspace(1,t-1,t-1),dk_wLC,'c');
hold on;
plot(linspace(1,t-1,t-1),dk_wLC2,'m');
title('k_t-k_{t-1}');
legend('LC','SVD','ML','wLC','wLC2');


subplot(6,2,3) 
autocorr(dk_LC)
title('Sample Autocorrelation Function for dk_{LC}')
subplot(6,2,4)
parcorr(dk_LC)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for dk_{LC}')

subplot(6,2,5) 
autocorr(dk_SVD)
title('Sample Autocorrelation Function for dk_{SVD}')
subplot(6,2,6)
parcorr(dk_SVD)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for dk_{SVD}')

subplot(6,2,7) 
autocorr(dk_ML)
title('Sample Autocorrelation Function for dk_{ML}')
subplot(6,2,8)
parcorr(dk_ML)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for dk_{ML}')

subplot(6,2,9) 
autocorr(dk_wLC)
title('Sample Autocorrelation Function for dk_{wLC}')
subplot(6,2,10)
parcorr(dk_wLC)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for dk_{wLC}')

subplot(6,2,11) 
autocorr(dk_wLC2)
title('Sample Autocorrelation Function for dk_{wLC2}')
subplot(6,2,12)
parcorr(dk_wLC2)
ylabel('Sample Part. Autocor.')
title('Sample Partial Autocorrelation Function for dk_{wLC2}')

clearvars dk_LC dk_SVD dk_ML dk_wLC dk_wLC2 dA
