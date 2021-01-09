% u1' = 3*u1 + 2*u2 - (2t^2 + 1)*e^2t , u1(0) = 1
% u2' = 4*u1 + u2+ (t^2 + 2t -4)*e^2t  , u2(0) = 1

t0 = 0; t1 = 1;  % 0 <= t <= 1 
w1_0 = 1;   % u1(0) = 1
w2_0 = 1;   % u2(0) = 1
h = 0.2;

% u1' = f1(t, u1, u2), u2' = f2(t, u1, u2)

syms f1(t,u1,u2) f2(t,u1,u2)
f1(t,u1,u2) = 3*u1 + 2*u2 - (2*t^2 + 1)*exp(2*t);
f2(t,u1,u2) = 4*u1 + u2 +(t^2 + 2*t -4)*exp(2*t);

[t,w] = RK4_system2(t0,t1,h,w1_0,w2_0,f1,f2);

disp("Approximate values")
disp(w)

range = t0:h:t1;
y1 = [];
y2 = []; 

% Generating values of the actual solution
for t = t0:h:t1 %iterate through t=0 to t=1 with step-size 0.2
    y1(end+1) = (1/3)*exp(5*t) - (1/3)*exp(-t)+ exp(2*t); %append values to y array
    y2(end+1) = (1/3)*exp(5*t) + (2/3)*exp(-t) + t^2*exp(2*t);
end

% Error between approximate and exact solutions - difference between row 1
% of w and y1, row 2 of w and y2 
error_1 = [];
error_2 = [];

%iterate over rows of w 
for i = 1:size(w,1)
    if i == 1
        for j = 1:size(w,2)   % iterate through cols 
            error_1(j) = abs(y1(j)-w(i,j)); 
        end
    else %row ==2
        for j = 1:size(w,2)
            error_2(j) = abs(y2(j)-w(i,j)); 
        end
    end  
end


disp("Exact solution:  u1(t)")
disp(y1)
disp("Absolute error: |u1(t_i) - w1(i)|")
disp(error_1)

disp("Exact solution:  u2(t)")
disp(y2)
disp("Absolute error: |u2(t_i) - w2(i)|")
disp(error_2)





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

u_1 = (1/3)*exp(5.*t) - (1/3)*exp(-1.*t)+ exp(2.*t);
u_2 = (1/3)*exp(5.*t) + (2/3)*exp(-1.*t) + t.*t.*exp(2.*t);

figure
%plot(t,w(1,:),'-o',t,u_1,'-*');
%legend('RK4','Exact u1(t)', 'Location', 'northwest');
plot(t,w(2,:),'-o',t,u_2,'-*');
legend('RK4','Exact u2(t)', 'Location', 'northwest');
end
