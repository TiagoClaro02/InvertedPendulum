%% ?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("\n### Exercise 1 ###\n");

%% Data

% System Variables
M_int = 5;  % kg
m_int = 1;  % kg
g_int = 10; % ms^-2
l_int = 1;  % m

% Symbolic Variables
syms M m g L real;

syms x x_d x_dd real;
syms theta theta_d theta_dd real;
syms x1 x2 x3 x4;
syms x1_d x2_d x3_d x4_d;
syms X_d;
syms y1 y2;
syms u;

% System Equations
x_dd = (m*L*(theta_d)^2*sin(theta) - m*g*sin(theta)*cos(theta) + u) / (M+m-m*cos(theta)^2);
theta_dd = ((M + m)*g*sin(theta) - m*L*(theta_d)^2*sin(theta)*cos(theta) - u*cos(theta)) / (L*(M + m - m*cos(theta)^2));

%% a)   %************************************************************* DONE 

disp("Part A");

% variables substituition (theta for x1 and theta_d for x2)
x_dd = subs(x_dd, [theta, theta_d], [x1, x2]);
theta_dd = subs(theta_dd, [theta, theta_d], [x1, x2]);

% Create state vector X'
%       |x1'|   x1' = x2
% X' =  |x2'|   x3' = x4
%       |x3'|   x2 -> theta'
%       |x4'|   x4 -> x'

x1_d = x2;
x2_d = theta_dd;
x3_d = x4;
x4_d = x_dd;
X_d = [x1_d; x2_d; x3_d; x4_d];
f = X_d;

% System output
%
% theta-> x1
% x-> x3

y1 = x1;
y2 = x3;
h = [y1; y2];

% Nonlinear state-space representation
display(f); % f(x,u)
display(h); % h(x,u)

%% b)   %************************************************************* DONE 

disp("Part B");

x1_eq = 0;
x2_eq = 0;
x3_eq = 0;
x4_eq = 0;

u_eq = 0;

y1_eq = x1_eq;
y2_eq = x3_eq;

% Substitute equilibrium values in the x'
f_eq = subs(X_d, [x1, x2, x3, x4, u], [x1_eq, x2_eq, x3_eq, x4_eq, u_eq]);

display(f_eq);

%% c)  %************************************************************* DONE

disp("Part C");

A = jacobian(f, [x1 , x2, x3, x4]);
B = jacobian(f, u);
C = jacobian(h, [x1, x2, x3, x4]);
D = jacobian(h, u);

A = subs(A, [x1, x2, x3, x4], [x1_eq, x2_eq, x3_eq, x4_eq]);
B = subs(B, x1, x1_eq);

display(A);
display(B);
display(C);
display(D);

%% ?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%

%   %************************************************************* DONE

fprintf("\n## Exercise 2 ###\n");

% Replacing variables with values
A = subs(A, [M, m, g, L], [M_int, m_int, g_int, l_int]);
B = subs(B, [M, m, g, L], [M_int, m_int, g_int, l_int]);

% Convert symbolic matrixes to double
A = double(A);
B = double(B);
C = double(C);
D = double(D);

display(A);
display(B);
display(C);
display(D);

% define system matrixes
sys = ss(A, B, C, D);

% Calculate eigenvalues
eigValues = eig(A);

display(eigValues);

% Check if system is stable, unstable or marginally stable
if (all(real(eigValues) < 0))
    disp("System is stable");
elseif (all(real(eigValues) <= 0))
    disp("System is asymptotic stable");
else
    disp("System is unstable");
end

%% ?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("\n## Exercise 3 ###\n");

%% a)   %************************************************************* DONE

disp("Part A");

% Create the controllability matrix
CtrbMatrix = ctrb(A, B);

display(CtrbMatrix);

% Calculate rank of controllability matrix
rankCtrbMatrix = rank(CtrbMatrix);
n = length(A);

display(rankCtrbMatrix);
display(n);

% Check if system is controllable
if (rankCtrbMatrix == n)
    disp("System is controllable");
else
    disp("System is not controllable");
end

%% b)   %************************************************************* DONE

fprintf("\nPart B\n");


% Define desired poles
P = [-2.00002 -2.00001 -1.9998 -1.9999];

% Calculate the gain matrix K
K = place(A, B, P);

display(K);

% Calculate the closed loop system matrixes
Acl = A - B*K;

% Calculate the eigenvalues of the closed loop system
Ecl = eig(Acl);

% Define the closed loop system
sys_cl = ss(Acl, B, C, D);

%% c)   %************************************************************* DONE

disp("Part C");

%% i)

disp("Part i");

x0 = [5*pi/180; 0.05; 0; 0];

% plot for 10 seconds
t = 0:0.01:10;

% Define the input
u = 0*t;

% Simulate the system
[y, t, ~] = lsim(sys_cl, u, t, x0);

% Plot the response
figure;
hold off;
subplot(1, 2, 1);
initial(sys_cl, x0, 10);
hold on;
title("Initial input x0 = [5*pi/180; 0.05; 0; 0]");
xlabel("Time (s)");
ylabel("Response (rad)");

%% ii)   %************************************************************* DONE

disp("Part ii");

x0 = [60*pi/180; 0.1; 0; 0];

% Simulate the system
[y, t, x] = lsim(sys_cl, u, t, x0);

% Plot the response
subplot(1, 2, 2);
initial(sys_cl, x0, 10);
hold on;
title("Initial input x0 = [60*pi/180; 0.1; 0; 0]");
xlabel("Time (s)");
ylabel("Response (rad)");

%% d)   %************************************************************* DONE

disp("Part D");

% Define the system parameters
dx_dt = @(t,x) [x(2);
    ((M_int + m_int)*g_int*sin(x(1)) - m_int*l_int*x(2)^2*sin(x(1))*cos(x(1)) - 0*cos(x(1))) / (l_int*(M_int + m_int - m_int*cos(x(1))^2));
    x(4);
    (m_int*l_int*x(2)^2*sin(x(1))- m_int*g_int*sin(x(1))*cos(x(1)) + 0) / (M_int + m_int - m_int*cos(x(1))^2 ) ;
    ];

time_interval = [0,10]; % time interval for simulation
x0 = [5*pi/180 0.05 0 0]'; % initial condition
[t, x] = ode45(dx_dt, time_interval, x0);

% Plot the response
figure;
subplot(1, 2, 1);
plot(t, x(:,1), "Color", "r", "DisplayName", "x0 = [5*pi/180; 0.05; 0; 0]", "LineWidth", 1.5);
hold on;
title("Nonlinear angle");
xlabel("Time (s)");
ylabel("Response (rad)");

subplot(1, 2, 2);
plot(t, x(:,3), "Color", "r", "DisplayName", "x0 = [5*pi/180; 0.05; 0; 0]", "LineWidth", 1.5);
hold on;
title("Nonlinear position");
xlabel("Time (s)");
ylabel("Response (m)");

x0 = [60*pi/180 0.1 0 0]'; % initial condition
[t, x] = ode45(dx_dt, time_interval, x0);

% Plot the response
subplot(1, 2, 1);
plot(t, x(:,1), "Color", "b", "DisplayName", "x0 = [60*pi/180; 0.1; 0; 0]", "LineWidth", 1.5);
hold on;

legend;

subplot(1, 2, 2);
plot(t, x(:,3), "Color", "b", "DisplayName", "x0 = [60*pi/180; 0.1; 0; 0]", "LineWidth", 1.5);
hold on;

legend;

%% e)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part E");

%% ?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("\n## Exercise 4 ###\n");

%% a)   %************************************************************* DONE

disp("Part A");

% Create the observability matrix
obsvMatrix = obsv(A, C);

display(obsvMatrix);

% Calculate rank of observability matrix
rankObsvMatrix = rank(obsvMatrix);

display(rankObsvMatrix);
display(n);

% Check if system is observable
if (rankObsvMatrix == n)
    disp("System is observable");
else
    disp("System is not observable");
end

%% b)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

fprintf("\nPart B\n");

%% i)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part i");

%% ii)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part ii");

%% iii)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part iii");

%% iv)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part iv");

%% Part 5

fprintf("\n## Exercise 5 ###\n");

%% a)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part A");

%% b)   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO

disp("Part B");

