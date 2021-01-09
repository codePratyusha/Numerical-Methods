t0 = 0; t1 = 1; % Define the interval
h = 0.5;       % Step size
w0 = 0; % Initial condition

syms f(t,y)
f(t,y) = t*exp(3*t) -(2*y);

[t,w] = RK4(t0,t1,h,w0,f);
y = (1/5)*t.*exp(3*t) - (1/25)*exp(3*t) + (1/25)*exp(-2*t)
figure
plot(t,w, '-o',t,y,'-x')
legend('RK4','Exact', "Location","bestoutside")


function [t, w] = RK4(t0,t1,h,w0,f)
t = t0:h:t1;
w = zeros(size(t));
w(1) = w0;
for i = 1:size(t,2)-1
    k_1 = h*f(t(i),w(i));
    k_2 = h*f(t(i)+0.5*h,w(i)+0.5*k_1);
    k_3 = h*f((t(i)+0.5*h),(w(i)+0.5*k_2));
    k_4 = h*f((t(i)+h),(w(i)+k_3));
    w(i+1) = w(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4);
end
end
