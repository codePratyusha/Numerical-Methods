 %% System of equations

%eqns

f1 = @(x) x(1)^3 + x(1)^2*x(2) - x(1)*x(3) + 6 ;
f2 = @(x) exp(x(1)) + exp(x(2)) - x(3);
f3 = @(x) x(2)^2 -2*x(1)*x(3) - 4;



%system of eqns
f= @(x) [x(1)^3 + x(1)^2*x(2) - x(1)*x3;...
          exp(x1) + exp(x(2)) - x(3); ...
          x(2)^2 -2*x(1)*x(3)-4];

% J is jacobian

J = @(x) [3*x(1)^2+2*x(1)*x(2)-x(3) x(1)^2 -x(1);...
          exp(x(1)) exp(x(2)) -1;...
          -2*(x(3)) 2*x(2) -2*x(1)];

% F(x) = [ f1 ; f2; f3] s
%F = {f1 ; f2;  f3 };  
F = [];


% g = f1(x)^2 + f2(x)^2 + f3(x)^2 
%g = @(x) f1(x).^2 + f2(x).^2 + f3(x).^2 ;
% grad_g = 2 * (J(x))^transpose * F(x) 
%grad_g = @(x) 2 * J' * F ; 

% Initial input  x0 = [1; 1; 1]



% calculate grad_g 

g = [];  % array of g to store g values per iter
z = []  %array of grad_g

x1 = [1; 1; 1]
g(1) =  (f1(x1))^2 + (f2(x1))^2+ (f3(x1))^2 

% calculate grad_g 
J(x1);
J(x1)';
tempF1 = f1(x1);
tempF2 = f2(x1);
tempF3 = f3(x1);
F = [tempF1; tempF2 ; tempF3];

grad_g = 2*J(x1)'*F;
%z[1] = norm(grad_g);


x_init = [1;1;1];
[final_x_sd, store_x_sd, store_grad_sd] = SteepestDescent(x_init, 1e-5, 100, f1, f2, f3, J);
[final_x_newt, store_x_newt, store_grad_newt] = Newton_Method_Systems(x_init, 1e-5, 100, f1, f2, f3, J);

format long
disp("Final x for Steepest Descent")
disp(final_x_sd);
disp("Final x for Newton's Method for Systems")
disp(final_x_newt);
%display(store_grad_g);

%Testing how well the solutions work : looking for x such that F(x) = 0 

tempF1_fin_newt = f1(final_x_newt);
tempF2_fin_newt = f2(final_x_newt);
tempF3_fin_newt = f3(final_x_newt);

tempF1_fin_sd = f1(final_x_sd);
tempF2_fin_sd = f2(final_x_sd);
tempF3_fin_sd = f3(final_x_sd);

F_fin_newt = [tempF1_fin_newt; tempF2_fin_newt; tempF3_fin_newt];
F_fin_sd = [tempF1_fin_sd; tempF2_fin_sd; tempF3_fin_sd];

disp("F_soln_newt");
disp(F_fin_newt);
disp("F_soln_sd");
disp(F_fin_sd);


plot(store_grad_sd);
hold on 
plot(store_grad_newt, '-*');
hold off
legend("norm ( grad\_g(x^k) )  values for Steepest Descent","norm ( grad\_g(x^k) ) for Newton" , "Location","northeast");
%% Steepest Descent algorithm for systems of equations
function [x, store_x, z0] = SteepestDescent(x0,tol,maxIter, f1, f2, f3, J)

k=1;            %initialize iteration counter
%eps=1;          %initialize error
%alpha1=0; %set iteration parameter
%lpha3=1;
Kmax = maxIter;
%x1 = [1; 1; 1];
x = x0;    %start with given input
store_x = {1,100};
g = []; %store g values per iter
z0 = [] ;%store grad_gk values per iter 

%Computation loop
while k<Kmax
    g(k) =  (f1(x))^2 + (f2(x))^2+ (f3(x))^2  ;
    g1 = g(k);
    tempF1 = f1(x);
    tempF2 = f2(x);
    tempF3 = f3(x);
    F_x = [tempF1; tempF2; tempF3];
    grad_gk = 2 * J(x)' * F_x; 
    z0(k) = norm(grad_gk);   %storing ||grad_gk||
    
    if(z0 == 0)
            display("Zero gradient!")
            return 
    end
   
   % Step 5 in book 
   
   z = grad_gk / norm(grad_gk);
   alpha1 = 0; 
   alpha3 = 1;
   temp3 = x - alpha3*z;
   g3 = f1(temp3)^2 + f2(temp3)^2 + f3(temp3)^2;
   
   %STEP 6, 7 & 8
   while(g3 >= g1)
       alpha3 = alpha3/2;
       temp3 = x - alpha3*z; 
       g3 = f1(temp3)^2 + f2(temp3)^2 + f3(temp3)^2;
       
        if(alpha3 < tol/2)
                display("No likely improvement...")
                return
        end      
   end
   
   %STEP 9 
   alpha2 = alpha3/2;
   temp2 = x - alpha2*z; 
   g2 = f1(temp2)^2 + f2(temp2)^2 + f3(temp2)^2;
   
   %STEP 10 
    h1 = (g2-g1)/alpha2;
    h2 = (g3-g2)/(alpha3-alpha2);
    h3 = (h2-h1)/alpha3;
   
    %STEP 11 
    alpha0 = 0.5*(alpha2-h1/h3);
    temp0 = x - alpha0*z;
    g0 = f1(temp0)^2 + f2(temp0)^2 + f3(temp0)^2;
    
    %STEP 12 
    if(g0 <= g3)
            a = alpha0;
            g = g0;
    else
            a = alpha3;
            g = g3;
    end
   
    %STEP 13
    x_out = x - a*z;
    x=x_out;
    store_x(1,k) = {x_out};    %storing x vectors
    
    %STORING 
    if(abs(g-g1) < tol)
            display("The procedure was successful!")
            return
    end
        
    k = k+1; 
    
end 
    display("Maximum number of iterations exceeded!")
    return
    
end   
   
   
   



 
%% Newton's Method for systems of equations

function[x, store_x, z0] = Newton_Method_Systems(x0,tol, maxIter, f1, f2, f3, J)
    
k = 1;
x = x0;
    
    %store for x-vectors along the way
    store_x = {1,100};
    
    while(k <= maxIter)
        display("iteration");
        disp(k);
        
        tempF1 = f1(x);
        tempF2 = f2(x);
        tempF3 = f3(x);
        F_x = [tempF1; tempF2; tempF3];
        J_x = J(x);
        grad_gk = 2 * J(x)' * F_x; 
        z0(k) = norm(grad_gk);   %storing ||grad_gk||
        
        
        
       
        
        %display ("F(x) = \n", fx)
        % display ("J(x) = \n", jx)
        % display("J(x) Inverse = \n", jx_inv)
        
        y = linsolve(J_x,-F_x);
        %y = y'
        
        x = x + y ;
        
        
        store_x(1,k) = {x};
        
        %disp( y)
        %disp(x)
        
        if(norm(y) < tol)
            disp("The procedure was successful!")
            return 
        end
        
        k = k + 1; 
    end
        
    disp("Max number of iterations surpassed. The procedure was unsuccessful!")
    return 
end  
   
   


