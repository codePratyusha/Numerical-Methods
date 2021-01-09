% y'' = y' + 2y + cos(x) , 0 <= x <= pi/2 
% y(0) = -0.3 , y(pi/2) = - 0.1 


% Converting the above equation into a system of equationss we get :
% A : u1' = u2 , u2' = u2 + 2*u1 + cos(x) , 0 <= x <= pi/2, u1(0) = -0.3 , u2(0) = 0 
% B : v1' = v2 , v2' = v2 + 2*v1 , 0 <= x <= pi/2, v1(0) = 0, v2(0) = 1 


%Solving for A 
t0 = 0; t1 = pi/2;  % 0 <= t <= pi/2 
u1_0 = -0.3;   % u1(0) = -0.3
u2_0 = 0;   % u2(0) = 0
h = pi/4;

% u1' = f1(t, u1, u2), u2' = f2(t, u1, u2)

syms f1_u(t,u1,u2) f2_u(t,u1,u2)
f1_u(t,u1,u2) = u2;
f2_u(t,u1,u2) = u2 + 2*u1 + cos(t);

% w_u holds approximations for y1(t) 
[t,w_u] = RK4_system2(t0,t1,h,u1_0,u2_0,f1_u,f2_u);

% Solving for B
v1_0 = 0;   % v1(0) = 0
v2_0 = 1;   % v2(0) = 1

% v1' = f1(t,v1, v2), u2' = f2(t, v1, v2)
syms f1_v(t,v1,v2) f2_v(t,v1,v2)
f1_v(t,v1,v2) = v2;
f2_v(t,v1,v2) = v2 + 2*v1; 

% w_v holds approximations for y2(t) 
[t,w_v] = RK4_system2(t0,t1,h,v1_0,v2_0,f1_v,f2_v);

disp(' First row: approximations for u1(t), second row: approximations for u2(t):')
disp(w_u)
disp(' First row: approximations for v1(t), second row: approximations for v2(t):')
disp(w_v)

% using w_i = w_u,i  + ((-0.1) - w_u,pi/2)/w_v,pi/2 * w_v,i , generate w_i,
% which are approximations of y(ti) 

w = [];

for i = 1:3
    w(i) = w_u(1,i) + (-0.1 - w_u(1,3))/w_v(1,3) * w_v(1,i)
end

format long
disp('w_i : approximations of y(t)')
disp(w)

y = [];

% Generating values of the actual solution
for t = t0:h:t1 %iterate through t=0 to t=2 with step-size h
    y(end+1) = (-1/10)*(sin(t) + 3*cos(t)); %append values to y array
end

% Error between approximate and exact solutions - difference between row 1
% of w and y1
error = [];

for i = 1:length(y)
     error(i) = abs(y(i)-w(i)); 
end


disp('Exact solution, y(ti)')
disp(y)
disp('Absolute error, | y(ti) - w_i |')
disp(error)

figure
x = t0:h:t1;
plot(x,w, '-o', x,y, '-*');
legend('RK4 and Linear Shooting with h= pi/4','Exact', 'Location', 'northwest');

% RK4 for system of 2 ODEs (or second order ODE)
function [t,w] = RK4_system2(t0,t1,h,w1_0,w2_0,f1,f2)
t = t0:h:t1;
w = zeros(2,size(t,2));
w(1,1) = w1_0;
w(2,1) = w2_0;

for i = 1:size(t,2)-1
    k11 = h*f1(t(i),w(1,i),w(2,i));
    k12 = h*f2(t(i),w(1,i),w(2,i));
    k21 = h*f1(t(i)+h/2,w(1,i)+k11/2,w(2,i)+k12/2);
    k22 = h*f2(t(i)+h/2,w(1,i)+k11/2,w(2,i)+k12/2);
    k31 = h*f1(t(i)+h/2,w(1,i)+k21/2,w(2,i)+k22/2);
    k32 = h*f2(t(i)+h/2,w(1,i)+k21/2,w(2,i)+k22/2);
    k41 = h*f1(t(i+1),w(1,i)+k31,w(2,i)+k32);
    k42 = h*f2(t(i+1),w(1,i)+k31,w(2,i)+k32);
    w(1,i+1) = w(1,i)+1/6*(k11+2*k21+2*k31+k41);
    w(2,i+1) = w(2,i)+1/6*(k12+2*k22+2*k32+k42);
end
end
