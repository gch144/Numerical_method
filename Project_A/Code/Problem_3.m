% Numerical Methods, project A No. 9
% Question 3
% This Script contains a general program for solving a system of n linear
% equations Ax = b using Gauss Seidel and Jacobi Iterative method

clc;clear; close all
%% System Description
% The System Matrix
A = [14,-1,-6,5;
     1,-8,-4,-1;
     1,-4,-12,-1;
     1,-1,-8,-16];
% The output vector
b = [10;0;-10;-20];  
% Initialize the solution
x = zeros(length(b),1);         

%%  Solve the system using Gauss Seidel and Jacobi method
[xj, Rj,Kj] = Gauss_Jacobi(A,b, x);
[xs, Rs,Ks] = Gauss_Seidel(A, b, x);

%% Compare the results of iterations plotting norm of the solution error
figure()
hold on
plot(Kj,Rj,'Linewidth',2)
plot(Ks,Rs,'Linewidth',2)
grid on
xlabel('Iterations (n)')
ylabel('Error')
title('Norm of the solution error')
legend('Gauss-Jacobi','Gauss-Seidel')

%% Solving System of Problem 2 using Gauss-Seidel and Jacobi methods
n = 10;     % Number of equations

% Solution of Problem 2(a) using Gauss_Seidel method
fprintf('\rSolution of Problem 2(a) using Gauss_Seidel method\n')
System_Matrix = 'a'; 
% Generate System Matrix A and Output vector b
[A,b] = System_AB(n,System_Matrix);
% Initialize the solution
x = zeros(length(b),1); 
[x, Rs,Ks] = Gauss_Seidel(A, b, x);
Solution_xa = x
fprintf('Solution_error = %g \n',Rs(end))
fprintf('iterations = %g \n',Ks(end))

% Solution of Problem 2(b) using Gauss_Seidel method
fprintf('\rSolution of Problem 2(b) using Gauss_Seidel method\n')
System_Matrix = 'b'; 
% Generate System Matrix A and Output vector b
[A,b] = System_AB(n,System_Matrix); 
% Initialize the solution
x = zeros(length(b),1); 
[x, Rs,Ks] = Gauss_Seidel(A, b, x);
Solution_xb = x
fprintf('Solution_error = %g \n',Rs(end))
fprintf('iterations = %g \n',Ks(end))