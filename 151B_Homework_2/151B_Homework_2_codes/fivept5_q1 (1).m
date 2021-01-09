a = 0; b = 1; alpha = 0;
hmax = 0.25; hmin = 0.05; 
TOL = 1e-4;
global T ;
T = [a];
global W ;
W = [alpha];

syms f(t,y)
f(t,y) = t*exp(3*t) -(2*y);

[t, w] = RKF(a, b, alpha, TOL, hmax, hmin, f);

disp(t)
disp(w)

exact = (1/5)*t.*exp(3*t) - (1/25)*exp(3*t) + (1/25)*exp(-2*t);
figure
plot(t,w, '-o',t,exact, '--x', 'linewidth', 1);
legend('RKF','Exact', "Location","bestoutside")
disp(['error between the approximation and the exact solution is: ', num2str(w-exact)])


function [T, W] = RKF(a, b, alpha, TOL, hmax, hmin, f)
t = a; w = alpha; h = hmax; FLAG = 1;
global T; global W;
while FLAG == 1
    k1 = double(h*f(t,w));
    k2 = double(h*f(t+1/4*h,w+1/4*k1));
    k3 = double(h*f(t+3/8*h,w+3/32*k1+9/32*k2));
    k4 = double(h*f(t+12/13*h,w+1932/2197*k1-7200/2197*k2+7296/2197*k3));
    k5 = double(h*f(t+h,w+439/216*k1-8*k2+3680/513*k3-845/4104*k4));
    k6 = double(h*f(t+1/2*h,w-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5));
    R = abs(1/h*(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6));
    if R < TOL
        t = t+h;
        w = w+25/216*k1+1408/2565*k3+2197/4104*k4-1/5*k5;
        T = [T,t];
        W = [W,w];
    end
    delta = 0.84*(TOL/R)^(0.25);
    if delta<0.1
        h = 0.1*h;
    else
        if delta>= 4
            h = 4*h;
        else
            h = delta*h;
        end
    end
    if h>hmax
        h = hmax;
    end
    if t>=b 
        FLAG = 0;
    else
        if t+h>b
            h = b-t;
        else
            if h<hmin
                % FLAG = 0;
                disp('minimum h exceeded')
                return
            end
        end
    end
    %disp('The procedure is complete')
end

end
