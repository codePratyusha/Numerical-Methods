%% Non-linear Shooting method: Secant method
% y'' = y^3 - y*y' , 1 <= x <= 2, y(1) = 1/2, y(2) = 1/3 
% h = 0.1
% tol = 1e-4

maxiter = 100;      
t = zeros(1,maxiter);   %store all our guesses of the initial derivative 
t(1) = -1/6 ; t(2) = -0.2415614935;   %set initial guesses : t(1) = theta_0 , t(2) = theta_1.  


t0 = 1; t1 = 2;     %boundary values
w0 = [0.5, t(1)];   % w0 is y(1)
h = 0.1;
tol = 1e-4;
yb = zeros(1, maxiter);  % storing shooting curves
beta = 1/3;

% y is a 1 x 2 vector : Let y(1) = y, y(2) = y' 
%Rewriting the given y'' gives
f = @(x,y) [y(2) ;  (y(1))^3-y(1)*y(2)];  % [y' ; y^3 - y*y']
[t_0, w_0] = RK4_system(t0, t1, h, w0, f); 

figure
disp("yb(t(1))")
disp(w_0(1,end))
plot(t_0, w_0(1,:)) 

hold on

yb(1) = w_0(1,end);
err = 100;     %initialize some error to start
i = 2;  %iterator

while(err > tol)% compute the shooting curve according to the initial guess t(i)
    w0 = [0.5, t(i)]; 
    [t_0, w_0] = RK4_system(t0, t1, h, w0, f);
    plot(t_0, w_0(1,:));
    yb(i) =  w_0(1,end);   %result of shooting curve
    t(i+1) = t(i) - (yb(i) - beta)*(t(i) - t(i-1))/(yb(i)-yb(i-1));  %update t(i+1)
    %update error
    disp("Errors, secant");
    err = abs(yb(i) - beta);
    disp(err);
    
    i = i+1;
end

% y_true is actual solution
y_true = 1./(t_0 + 1);
plot(t_0,y_true, '-*');
hold off


% Error between approximate and exact solutions - difference between row 1
% of w and y1
error_sec = [];
t_0_sec = t_0; %copy x values to store
for i = 1:length(y_true)
     error_sec(i) = norm(w_0(1,i)-y_true(i)); 
end

format long
% disp('For iter = 10, the approximations for y(t_i), w_i , is :')
% disp(w_0(1,:))
% disp('Exact solution y(ti):')
% disp(y_true)
% disp('Absolute error, | y(ti) - w_i |')
% disp(error)

%iteratively add info to legend
cell(1, 3);
for i = 1:3
    str{i} = sprintf('t_%d = %0.7f' , i , t(i));
end
str{4} = 'Exact';
legend(str, 'Location', 'southoutside')

%Secant method : t_(k+1) = t_k - f(f(k-beta)*(t_k-t_(k-1))/(f(k)-f(k-1))
% we want f(k) to approach beta, the boundary value of y(t1)



%% Non-linear Shooting method: Newton's method

% y'' = y^3 - y*y' , 1 <= x <= 2, y(1) = 1/2, y(2) = 1/3 
% h = 0.1
% tol = 1e-4

tol = 1e-4;
maxiter = 100;
t = zeros(1,maxiter); % store all our guesses of the initial derivatives
t(1) = -1/6;          % initial guesses, t(1) = theta_0

t0 = 1; t1 = 2; 
w0 = [0.5, t(1)];     % [y(1); t(1)]
h = 0.1;

yb = zeros(1,maxiter); 
beta = 1/3;

% y is a 1 by 2 vector: y(1) = y, y(2) = y'
% z is a 1 by 2 vector: z(1) = z, z(2) = z'
f = @(x,y) [y(2); (y(1))^3-y(1)*y(2)];
g = @(x,z,y) [z(2); z(1)*(3*y(1)^2 - y(2))-y(1)*z(2)]; % computes the derivative of z
[t_0,w_0] = RK4_system(t0,t1,h,w0,f);
figure
plot(t_0,w_0(1,:))
hold on
yb(1) = w_0(1,end);
err = 100; %set arbitrary value for err
i = 2;   %iterator

while(err > tol) % compute the shooting curve according to the initial guess t(i)
    w0 = [0.5, t(i)]; z0 = [0,1];
    [t_0,w_0] = RK4_system(t0,t1,h,w0,f);
    plot(t_0,w_0(1,:))
    yb(i) = w_0(1,end);
    [t_0,z_0] = RK4_newton(t0,t1,h,z0,g,w_0);
    %Newton's method: t(i+1) = t(i)-(yb(i)-beta)/(dydt at b) (:equals z_0(end))
    t(i+1) = t(i)-(yb(i)-beta)/z_0(1,end);    %update t(i+1)
    
    %update error
    %disp("Errors");
    err = abs(yb(i) - beta);
    %disp(err);
    
    i = i+1;
end

% y_true is actual solution
y_true = 1./(t_0 + 1);
plot(t_0,y_true, '-*');
hold off

% Error between approximate and exact solutions - difference between row 1
% of w and y_true
error_newt = [];   % holds error
t_0_newt = t_0;
for i = 1:length(y_true)
     error_newt(i) = norm(w_0(1,i)- y_true(i)); 
end


str = cell(1,maxiter);

for i = 1:maxiter
    str{i} = sprintf('t_%d = %0.7f',i,t(i)) ;
end
str{5} = "Exact"; 
legend(str, "Location",'southoutside')


%% Error plots 

figure
hold on 
semilogy(t_0_sec, error_sec, '-o');
semilogy(t_0_newt, error_newt, '-*');
legend("secant err", "newton err", 'Location','southoutside')
hold off


%% RK4 for newton shooting
function [t,z] = RK4_newton(t0,t1,h,z0,f,y)
t = t0:h:t1;
z = zeros(length(z0),length(t)); % initialize w as length(t) vectors. each vector stores w1 and w2 for 2-system

z(:,1) = z0;
for i = 1:length(t)-1
    k1 = h*f(t(i),z(:,i),y(:,i));
    k2 = h*f(t(i)+h/2,z(:,i)+k1/2,y(:,i));
    k3 = h*f(t(i)+h/2,z(:,i)+k2/2,y(:,i));
    k4 = h*f(t(i+1),z(:,i)+k3,y(:,i));
    z(:,i+1) = z(:,i)+1/6*(k1+2*k2+2*k3+k4);
end
end

%% RK4 for system of ODEs
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




