%% Q4 helper code, mostly done by hand 


A = [-1 -2 1 2; -2 3 0 -2; 1 0 2 1; 2 -2 1 4];

a_21 = -2;
a_31 = 1;
a_41 = 2;

alpha = -1*sqrt(a_21^2 + a_31^2 + a_41^2);
r = sqrt(1/2*(alpha)^2 - 1/2*(a_21)*(alpha));

w1 = 0;
w2 = (a_21 - alpha)/(2*r);
w3 = (a_31)/(2*r);
w4 = (a_41)/(2*r);

w = [w1; w2; w3; w4];

P1 = eye(4) - 2*w*w';

A2 = P1*A*P1;

disp("alpha")
disp(alpha);
disp("r");
disp(r);
disp("w");
disp(w);
disp("P1");
disp(P1);
disp("A2");
disp(A2);

a_32 = 0.8889;
a_42 = 0.1111;

alpha_2 = -1*sqrt(a_32^2 + a_42^2);
disp("alpha_2")
disp(alpha_2);

r_2 = sqrt(1/2*(alpha_2)^2 - 1/2*(a_32)*(alpha_2));
disp("r_2");
disp(r_2);

w1 = 0; 
w2 = 0 ;
w3 = (a_32 - alpha_2)/(2*r_2);
w4 = a_42/(2*r_2);

w = [w1; w2; w3; w4];
disp("w");
disp(w);


P2 = eye(4) - 2*w*w';
disp("P2");
disp(P2);


A3 = P2*A2*P2;
disp("A3");
disp(A3);