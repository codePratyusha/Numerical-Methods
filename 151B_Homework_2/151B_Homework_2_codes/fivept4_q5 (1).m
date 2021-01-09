% Midpoint: k1 = h/2*f(t(i),w(i)) : approximate the derivative of the
% midpoint, which is (t(i)+t(i+1))/2
% w(i+1) = w(i)+h*f(t(i)+h/2,w(i)+k1/2) : w(i)+k1/2 is the approximated
% value of y(t(i)+h/2) (which is the midpoint)

t0 = 0; t1 = 1; % Define the interval
h = 0.5;       % Step size
w0 = 0; % Initial condition

syms f(t,y)
f(t,y) = t*exp(3*t) -(2*y);

[t,w] = Midpoint(t0,t1,h,w0,f);
y = (1/5)*t.*exp(3*t) - (1/25)*exp(3*t) + (1/25)*exp(-2*t)
figure
plot(t,w, '-o',t,y,'-x')
legend('Midpoint','Exact', "Location","bestoutside")

 error = zeros(1,5);
 yexact = (1/5)*t1.*exp(3*t1) - (1/25)*exp(3*t1) + (1/25)*exp(-2*t1) 
 for i = 1:5
    [t,w] = Midpoint(t0,t1,h,w0,f);
    error(i) = abs(w(end) - yexact)
end    

loglog(error(1:4))
xlabel('step size')
ylabel('error')


function [t,w] = Midpoint(t0,t1,h,w0,f)
t = t0:h:t1;
w = zeros(size(t));
w(1) = w0;
for i = 1:size(t,2)-1
    k1 = h*f(t(i),w(i));
    w(i+1) = w(i)+h*f(t(i)+h/2,w(i)+k1/2);
end
end
