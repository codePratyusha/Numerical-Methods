% k1 = f(t(i+1),w(i)+h*f(t(i),w(i));
% w(i+1) = w(i)+h*(f(t(i),w(i))+k1)/2

t0 = 0; t1 = 1; % Define the interval
h = 0.5;       % Step size
w0 = 0;% Initial condition
syms f(t,y)
f(t,y) = t*exp(3*t) -(2*y);

[t,w_3] = Mod_Euler(t0,t1,h,w0,f);

y = (1/5)*t.*exp(3*t) - (1/25)*exp(3*t) + (1/25)*exp(-2*t)
figure
plot(t,w_3, '-o',t,y,'-x')
legend('Mod\_Euler','Exact', "Location","bestoutside")

error3 = zeros(1,5);
yexact = (1/5)*t1.*exp(3*t1) - (1/25)*exp(3*t1) + (1/25)*exp(-2*t1) 
for i = 1:5
     error3(i) = abs(w_3(end) - yexact)
end    

loglog(error3(1:4))
xlabel('step size')
ylabel('error')



function [t,w] = Mod_Euler(t0,t1,h,w0,f)
t = t0:h:t1;
w = zeros(size(t));
w(1) = w0;
for i = 1:size(t,2)-1
    k1 = f(t(i+1),w(i)+h*f(t(i),w(i)));
	w(i+1) = w(i)+h*(f(t(i),w(i))+k1)/2;
end
end





