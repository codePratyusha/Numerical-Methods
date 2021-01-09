% Non-linear Shooting method: Secant method for Q11.2.1
% y'' = -(y')^2 -y + ln(x)  , 1<= x <= 2, y(1) = 0, y(2) = ln(2) 
%h  = 0.5
maxiter = 10 ;
t = zeros(1,maxiter);   %store all our guesses of the initial derivative 
t(1) = 0.2 ; t(2) = 0.4;   %arbitrarily change initial guesses

t0 = 1; t1 = 2; 
w0 = [0, t(1)];   %
h = 0.5;
yb = zeros(1, maxiter);
beta = log(2);

% y is a 1 x 2 vector : y(1) = y, y(2) = y'
f = @(x,y) [y(2) ; -(y(2))^2 - y(1) + log(x)];  % y(1) = y, y(2) = y'
[t_0, w_0] = RK4_system(t0, t1, h, w0, f);


figure
%disp(w_0(1,end))
plot(t_0, w_0(1,:))

hold on

yb(1) = w_0(1,end);
for i = 2:maxiter % compute the shooting curve according to the initial guess t(i)
    w0 = [0, t(i)]; 
    [t_0, w_0] = RK4_system(t0, t1, h, w0, f);
    plot(t_0, w_0(1,:));
    yb(i) =  w_0(1,end);   %result of shooting curve
    t(i+1) = t(i) - (yb(i) - beta)*(t(i) - t(i-1))/(yb(i)-yb(i-1));  %update t(i+1)
end
y = log(t_0);
plot(t_0,y, '-*');

hold off



% Error between approximate and exact solutions - difference between row 1
% of w and y1
error = [];

for i = 1:length(y)
     error(i) = abs(y(i)-w_0(1,i)); 
end

format long
disp('For iter = 10, the approximations for y(t_i), w_i , is :')
disp(w_0(1,:))
disp('Exact solution y(ti):')
disp(y)
disp('Absolute error, | y(ti) - w_i |')
disp(error)

%iteratively add info to legend
cell(1, maxiter);
for i = 1:maxiter
    str{i} = sprintf('t_%d = %0.7f' , i , t(i));
end
tr{maxiter+1} = 'Exact';
legend(str, 'Location', 'northwest');
%Secant method : t_(k+1) = t_k - f(f(k-beta)*(t_k-t_(k-1))/(f(k)-f(k-1))
% we want f(k) to approach beta, the boundary value of y(t1)


% RK4 for system of ODEs
function [t,w] = RK4_system(t0,t1,h,w0,f)
t = t0:h:t1;
w = zeros(length(w0),length(t)); % initialize w as length(t) vectors. each vector stores w1 and w2 for 2-system

w(:,1) = w0;
for i = 1:length(t)-1
    k1 = h*f(t(i),w(:,i));
    k2 = h*f(t(i)+h/2,w(:,i)+k1/2);
    k3 = h*f(t(i)+h/2,w(:,i)+k2/2);
    k4 = h*f(t(i+1),w(:,i)+k3);
    w(:,i+1) = w(:,i)+1/6*(k1+2*k2+2*k3+k4);
end
end

