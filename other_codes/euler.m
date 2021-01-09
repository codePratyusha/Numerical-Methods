t0 = 0; t1 = 1; % Define the interval
h = 0.1;       % Step size
w0 = sqrt(2);         % Initial condition

syms f(t,y)
f(t,y) = (50/y) - (50*y);
w = Euler1(t0,t1,h,w0,f);

function [t,w] = Euler1(t0,t1,h,w0,f)
t = t0:h:t1;
w = zeros(size(t));
w(1) = w0;
for i = 1:size(t,2)-1
    w(i+1)=w(i)+h*f(t(i),w(i));
end

y = sqrt(1 + exp(-100*t));
figure
plot(t,w,'-o',t,y,'-*')
legend('Euler','Exact', 'Location','northwest')

end
