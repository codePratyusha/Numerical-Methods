% u1' = u2 , 
% u2' = (2/t)*u2 - (2/t^2)*u1 + t*log(t) 
% u1(1) = 1 , u1'(0) = 0 = u2(1)

t0 = 1; t1 = 2;  % 1 <= t <= 2 
w1_0 = 1;   % u1(1) = 1
w2_0 = 0;   % u2(1) = 0
h = 0.1;

% u1' = f1(t, u1, u2), u2' = f2(t, u1, u2)

syms f1(t,u1,u2) f2(t,u1,u2)
f1(t,u1,u2) = u2;
f2(t,u1,u2) = (2/t)*u2 - (2/t^2)*u1 + t*log(t) ;

[t,w] = RK4_system2(t0,t1,h,w1_0,w2_0,f1,f2);

format long
disp("Approximate values")
disp(w.')

y = [];

% Generating values of the actual solution
for t = t0:h:t1 %iterate through t=0 to t=2 with step-size h
    y(end+1) = (7/4)*t + (1/2)*t^3*log(t) - (3/4)*t^3; %append values to y array
end

% Error between approximate and exact solutions - difference between row 1
% of w and y1
error_1 = [];

%iterate over rows of w 
for i = 1:size(w,1)
    if i == 1
        for j = 1:size(w,2)   % iterate through cols 
            error_1(j) = abs(y(j)-w(i,j)); 
        end
    end  
end



disp("Exact solution:  y(t)")
disp(y.')
disp("Absolute error: |y(t_i) - w1(i)|")
disp(error_1.')


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

y = (7/4).*t + (1/2).*t.*t.*t.*log(t) - (3/4).*t.*t.*t ;

figure
plot(t,w(1,:),'-o',t,y,'-*');
 legend('RK4','Exact', 'Location', 'southwest');

end
