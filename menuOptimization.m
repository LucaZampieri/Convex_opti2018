clear all
close all
clc

%% Set the parameters for the barrier method
t = 1;
mu = 20;
tol = 1e-10;
maxIter = 100;

tolNewton = 1e-5;
maxIterNewton = 100;

%% Define the objective function and the constraints
A = [199, 11.37, 6.71, 16.72; 191, 13, 5.2, 13; 128, 5.5, 5.5, 14; 191, 7.43, 17.36, 12.58; 116, 1.88, 21.04, 4.14]; % Nutrient constraints
C = [5.5, 4.5, 5.0, 3.5, 3.4]'; % Cost function
A = [C, A]';
[m, n] = size(A);
A = [A; -A; ones(1, n); -ones(1, n); -eye(n)];
b2 = [14.0, 450, 30, 50, 60]';
b1 = [0, 0, 0, 0, 30]';
food_max = 3;
food_min = 0.8;
b = [b2; -b1; food_max; -food_min; zeros(n,1)];

T = [4.9, 3.9, 4.4, 4, 4.4]'; %Taste function
x = 0.5*ones(n, 1);

hess_f0 = @(x) x * zeros(1, n);
phi = @(x) -sum(log(b-A*x));
grad_phi = @(x) A' * (1./(b-A*x));
hess_phi = @(x) A' * diag((1./(b-A*x)).^2) * A;

%% Solve the otpimization problem for different values of gamma
gammaMin = 0.0;
gammaMax = 1.6;
N = 10;
gamma = linspace(gammaMin, gammaMax, N);
xOpt = zeros(n, N);
for i = 1:N
    
    f0 = @(x) (gamma(i) * C' - T') * x;
    grad_f0 = @(x) gamma(i) * C -  T;

    xOpt(:, i) = barrierMethod(f0, grad_f0, hess_f0, phi, grad_phi, hess_phi, x, t, mu, m, tol, maxIter, tolNewton, maxIterNewton);
    i
end

    
    
    
    
