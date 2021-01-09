% RK4 for first order ODE.

t0 = 0; t1 = 1; % Define the interval
h = 0.1;       % Step size
w0 = sqrt(2);         % Initial condition

syms f(t,y)
f(t,y) = (50/y) - (50*y);
w = RK4(t0,t1,h,w0,f);

function [t,w] = RK4(t0,t1,h,w0,f)
t = t0:h:t1;
w = zeros(size(t));
w(1) = w0;
for i = 1:size(t,2)-1
    k1 = h*f(t(i),w(i));
    k2 = h*f(t(i)+h/2,w(i)+k1/2);
    k3 = h*f(t(i)+h/2,w(i)+k2/2);
    k4 = h*f(t(i+1),w(i)+k3);
    w(i+1) = w(i)+1/6*(k1+2*k2+2*k3+k4);
end
y = sqrt(1 + exp(-100*t));
figure
plot(t,w,'-o',t,y,'-*')
legend('RK4','Exact', 'Location','southwest')
end
