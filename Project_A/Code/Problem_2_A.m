% Numerical Methods, project A No. 9
% Question 2
% This Script contains a general program for solving a system of n linear
% equations Ax = b using Gaussian elimination with partial pivoting

% Clear and clean up the works space
clc;clear;close all
format shortg
format compact
% Set the value of n
n0 = 10; % Initial number of equations

%% Problem 2 case (a)
disp('Problem 2 case (a)')
System_Matrix = 'a'; 
% Iterate for increasing values of
% n = 10,20,40,80,160,… until the solution time becomes prohibitive
 n(1) = n0; % Start at n=10
 ErrorNorm_a = zeros(1,20); % Initialize the Error norm
for i=1:4%8  
    % Generate the system matrix and output Vector b
    [A,b] = System_AB(n(i),System_Matrix);
    x = GaussPP(A,b);
    % Calculate the Euclidean norm of residuum r = Ax–b
    ErrorNorm_a(i) = norm(A*x - b);
    n(i+1) = 2*n(i);    
end
ErrorNorm_a = ErrorNorm_a(1:i);
% plot Error versus n.
plot(n(1:end-1),ErrorNorm_a)
grid on
xlabel('Number of equation [n]')
ylabel('Error norm')
title('Solution error vs n')
% For n = 10 print the solutions and the solutions’ errors
fprintf('Solutions and the solutions’ errors for n = 10\n')
[A,b] = System_AB(10,'a'); % Get matrix A and vector b
Solution_x = GaussPP(A,b);           % System Solution
solutions_errors = A*Solution_x - b; % Solution Errors
T = table(Solution_x,solutions_errors)

% Making the residual correction and checking if it improves the solutions.
fprintf('\rSolutions and the solutions’ errors for n = 10 after residual correction\n')
err = GaussPP(A,solutions_errors);
Solution_x = Solution_x-err;
solutions_errors = A*Solution_x - b;
T = table(Solution_x,solutions_errors)
