%% Numerical Methods, PROJECT C No. 10
% Problem I:
% We are required to determine a polynomial function y=f(x) that best fits
% the experimental data by using the least-squares approximation by testing
% polynomials of various degrees

% Clear and prepare the Matlab workspace
clc;clear all;close all
% The following are the experimental measurements (samples):
x = (-5:5);
y = [-4.6643,-5.5445,-4.0378,-2.0868,-1.0393,0.6916,-0.0237,0.4107,...
    -3.6761,-9.8466,-18.8868];


%% Least Square approximation
% To solve the least-squares problem use the system of normal equations
% with QR factorization of a matrix A.For each solution calculate the error
% defined as the Euclidean norm of the vector of residuum and the condition
% number of the Gram's matrix. Compare the results in terms of solutions’
% errors

N = length(x)-1; % Highest polynomial degree depends of number of samples
y_fit = zeros(length(x),N); % Initialize variable to store fitted values
CNum = zeros(N,1); % Initialze a varible to store condition numbers
% Euclidean norm of the vector of residuum with A matrix
ENorm_A = zeros(N,1);
% Euclidean norm of the vector of residuum with Gram's Matrix

ENorm_G = zeros(N,1);
for degree_n=1:N
    % Generate a system of normal equations for polynomial order j
    [A,b]= SystemOfEquations(x,y,degree_n);
    % solve the least-squares problem use the system of normal equations
    % with QR factorization of a matrix A
    LS_Sol = QR_solution(A,b);% Least Square (LS) Solutions
    y_fit(:,degree_n) = polyval(fliplr(LS_Sol(:)'),x(:));
    
    % Euclidean norm of the vector of residuum with Gram's Matrix
    GM = A'*A;                              % Gram's matrix
    CNum(degree_n) = cond(GM);              % Cond No.of Gram's matrix
    residuum_gm = GM*LS_Sol(:) - (A'*b);    % Residuum
    ENorm_G(degree_n) = norm(residuum_gm,2);% Euclidean norm
    
    % Euclidean norm of the vector of residuum with A matrix
    residuum_A = A*LS_Sol(:) - b;           % Residuum
    ENorm_A(degree_n) = norm(residuum_A,2); % Euclidean norm
    clear coeff A
end
% Display the results in a table format
PolynomialDegree=[1:10]';
Results = table(PolynomialDegree,ENorm_A,CNum, ENorm_G,...
    'VariableNames',...
    {'Polynomial Degree', 'Euclidean Norm of residuum of A matrix',...
    'Condition number of Gram''s matrix', ...
    'Euclidean Norm of residuum of Gram''s Matrix'})

%% Plot the results
% Plot the solutions for of the polynomials of various degrees)
p_order ={'First','Second','Third','Fourth','Fifth','Sixth','Seventh',...
    'Eigth','Ninth','Tenth'};
for j=1:N
    figure(); % Generate a new figure window for each polynomial degree
    hold on % Use hold to place different plots on the same axis
    plot(x,y,'r*','markersize',10); % Plot the data samples
    plot(x,y_fit(:,j),'-b','linewidth',2)% Plot the fitted data
    xlabel('x');ylabel('y')% Add axis lables 
    grid on; % Add grid lines to the plot    
    % Add a title heading to the plot
    TitleHead = ['Least Square fit : ',p_order{j},' order polynomial'];
    title(TitleHead)
    % Add legend for each plotted curve
    legend('data samples','Least Square fit')
    hold off % Allos the the plotted curves to shown in the figure window
end

% Plot the condition number associated with the matrix of each polynomial
% order
figure()                                    % Create a new figure window
plot(1:N,CNum,'linewidth',2)                % Plot the condition number
xlabel('polynomial order ')                 % Add x-axis label
ylabel('Condition number')                  % Add y-axis label
grid on                                     % Add grid lines to the plot
title('Condition number of Gram''s matrix') % Add title of graph
 
% Plot Euclidean Norm of residuum of A matrix
figure()
plot(1:N,ENorm_A,'linewidth',2)                 % Add the plot
xlabel('polynomial degree ')                    % Add x-axis label
ylabel('Error')                                 % Add y-axis label
grid on                                         % Add grid lines on plot
title('Euclidean Norm of residuum of A matrix') % Add title of graph

% Plot the Error: Euclidean Norm of residuum of Gram's matrix
figure()                                        % Create a new figure window
plot(1:N,ENorm_G,'linewidth',2)                 % Add the plot
xlabel('polynomial degree ')                    % Add x-axis label
ylabel('Error')                                 % Add y-axis label
grid on                                         % Add grid lines on plot
title('Euclidean Norm of residuum of Gram''s matrix')% Add title

%% Functions
function [A,b]=SystemOfEquations(x,y,Poly_degree)
x = x(:); % data samples
b = y(:); % function values at data sample points
% Initialize the coefficients matrix of normal equations
A = [ones(size(x))]; 
for i=1:Poly_degree
    A(:,i+1)=x.^i;
end
end
function coeff=QR_solution(A,b)
b = b(:);        % Ensure it a a column vector
[~,n] = size(A); % Get the number of columns of matrix A
[Q,R,P] = qr(A); % Apply QR factorization
d = Q'*b;        % Compute matrix-vector product

% Apply back substitution
x = zeros(n,1); % Initialize vector x
for i = n:-1:1    
    x(i) = d(i);
    for j = i+1:n
        x(i) = x(i) - R(i,j)*x(j);
    end
    x(i) = x(i)/R(i,i);
end
coeff = P*x; % Re-order the coefficients using permutation matrix from QR
end
