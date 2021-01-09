
% [t0,t1] = interval of h 
% h = step size
% w0 = initial condition
% f = f(t,y) = y' 

t0 = 0; t1 = 1; % Define the interval
h = 0.1;       % Step size
w0 = sqrt(2);         % Initial condition

syms f(t,y)
f(t,y) = (50/y) - (50*y);

w = ABM4_PC(t0,t1,h,w0,f);

function [t,w] = ABM4_PC(t0,t1,h, w0,f)
% initialize arrays for plots
t = t0:h:t1;
w = zeros(size(t));

%initialize values
t(1) = t0;
w(1) = w0;

%generate values using RK4 
for i = 1:3
    k1 = h*f(t(i),w(i));
    k2 = h*f(t(i)+h/2,w(i)+k1/2);
    k3 = h*f(t(i)+h/2,w(i)+k2/2);
    k4 = h*f(t(i+1),w(i)+k3);
    w(i+1) = w(i)+1/6*(k1+2*k2+2*k3+k4);
end

% iterate
for i = 4:size(t,2)-1
	% Predict using AB
	w(i+1) = w(i) + h/24*(55*f(t(i), w(i)) - 59*f(t(i-1),w(i-1)) ...
		+ 37*f(t(i-2),w(i-2)) - 9*f(t(i-3),w(i-3)));
	t(i+1) = t(i) + h;
	% Correct using AM 
	w(i+1) = t(i) + h/24*(9*f(t(i+1),w(i+1)) + 19*f(t(i),w(i)) ...
		- 5*f(t(i-1),w(i-1)) + f(t(i-2),w(i-2)));
end
y = sqrt(1 + exp(-100*t));
figure
plot(t,w,'-o',t,y,'-*')
legend('ABM4_PC','Exact', 'Location','southwest')
end




