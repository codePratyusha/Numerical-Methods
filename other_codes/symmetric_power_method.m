%% a) Dominant eigenvalue

% Compute Eigenvalues in Matlab
A = [4 1 -1 0; 1 3 -1 0; -1 -1 5 2; 0 0 2 4];
display(A)
[V, D] = eig(A);  % diagonal entries of D are e-values, column vectors of V are e-vectors s.t. A*V = V*D
display(D)

%% Power Method with scaling : compute the largest eigenvalue and its corresponding eigenvector
w = zeros(4,4); % matrix of eigenvector approx. - kth col = kth iter e-vec approx --- 0: w(1); w(2); ...
x_approx = zeros(4,4);  % matrix of normalized eigenvector approx. -- same size as A -- x(0); x(1); ...
mu = []; % array of e-vals, 0, mu(1), mu(2), ... where mu(1) is dominant e-val


A = [4 1 -1 0; 1 3 -1 0; -1 -1 5 2; 0 0 2 4];
x0 = [0; 1; 0; 0];          %initialize with random non-zero col vector 
x0 = x0/norm(x0,"inf");           %normalize initalizer

%compute first iter with given x0
w(:,1) = A*x0 ; 
x_approx(:,1) = w(:,1)/ norm(w(:,1), "inf") ;
mu(1) = norm(w(:,1), "inf")

% disp("Norm w1")
% disp(norm(w( :,1)));
% disp(w( :,1))

%compute second iter 
w(:,2) = A*x_approx(:,1) ;             
x_approx(:,2) = w(:,2)/ norm(w(:,2), "inf") ;
mu(2) = norm(w(:,2), "inf")

% disp("Norm w2")
% disp(norm(w(:,2)));
% disp(w( :,2))

iter = 2;                    %set iter counter- iter updates once it has finished that iteration
k= 3;    
                                                                                      
while (norm(mu(iter) - mu(iter-1),1) > 10^(-4)) || (norm(x_approx(:, iter) - x_approx(:, iter-1), "inf") > 10^(-4)) || (iter <= 25 )
    
    w(:,k) = A*x_approx(:,k-1);                 %eigenvector approx
    x_approx(:,k) = w(:,k)/ norm(w(:,k), "inf");       % normalize e-vec approx 
    mu(k) = norm(w(:,k), "inf")     % e-val approx  using Rayleigh Quotient

    k = k+1 ;
    iter = iter + 1;
    
end


%display ("Eigenvectors approx")
%display(x_approx)
%display("Iterations")
%display(iter)
%plot(mu)
%display (norm(x_approx(:,2) - x_approx(:,1)));
%display(norm(mu(iter) - mu(iter-1)));


%% Symmetric Power Method : compute the largest eigenvalue and its corresponding eigenvector
w_sym = zeros(4,4); % matrix of eigenvector approx. - kth col = kth iter e-vec approx --- 0: w(1); w(2); ...
x_approx_sym = zeros(4,4);  % matrix of normalized eigenvector approx. -- same size as A -- x(0); x(1); ...
mu_sym = []; % array of e-vals, 0, mu(1), mu(2), ... where mu(1) is dominant e-val


A = [4 1 -1 0; 1 3 -1 0; -1 -1 5 2; 0 0 2 4];
x0 = [0; 1; 0; 0];          %initialize with random non-zero col vector 
x0 = x0/norm(x0, 2);           %normalize initalizer - 2-norm

%compute first iter with given x0
w_sym(:,1) = A*x0 ;             
x_approx_sym(:,1) = w_sym(:,1)/ norm(w_sym(:,1), 2) ;
mu_sym(1) = (x_approx_sym(:,1)'*A*x_approx_sym(:,1))

%compute second iter  
w_sym(:,2) = A*x_approx_sym(:,1) ;             
x_approx_sym(:,2) = w_sym(:,2)/ norm(w_sym(:,2), 2) ;
mu_sym(2) = (x_approx_sym(:,2)'*A*x_approx_sym(:,2))

iter = 2;                    %set iter counter- iter updates once it has finished that iteration
k= 3;    
                                                                                      
while (norm(mu_sym(iter) - mu_sym(iter-1),2) > 10^(-4)) || (norm(x_approx_sym(:, iter) - x_approx_sym(:, iter-1),2) > 10^(-4)) || (iter <= 25 )
    
    w_sym(:,k) = A*x_approx_sym(:,k-1);                 %eigenvector approx
    x_approx_sym(:,k) = w_sym(:,k)/ norm(w_sym(:,k),2);       % normalize e-vec approx 
    mu_sym(k) = (x_approx_sym(:,k)'*A*x_approx_sym(:,k))    % e-val approx  using Rayleigh Quotient

    k = k+1 ;
    iter = iter + 1;
    
end


format long
display("Eigenvalues approx for symmetric power method") 
display(mu_sym')
display("Eigenvalues approx for power method") 
display(mu')
%display ("Eigenvectors approx for sym")
%display(x_approx_sym)
%display("Iterations for sym ")
%display(iter)
%display (norm(x_approx(:,2) - x_approx(:,1)));
%display(norm(mu(iter) - mu(iter-1)));
plot(mu,'-o');
hold on
plot(mu_sym, '-*');
hold off
legend("Power Method" , "Symmetric Power Method", "Location", "southeast");
