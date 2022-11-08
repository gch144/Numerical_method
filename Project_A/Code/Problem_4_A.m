%This Scripts uses  QR method for finding eigenvalues of 5×5 matrices
%without shifts.
clc;clear;close all

% Chosen symmetric matrix 5×5
A = [10     9    13    15    16
     9    10     9    13    15
    13     9    10     9    13
    15    13     9    10     9
    16    15    13     9    10];

[m,n] = size(A);% m = number of Rows,  n = number of columns

% Initialize various variables with zeros
Q = zeros(m,n);
R = zeros(n,n);
L = tril(A,-1);
Lv = L(:); B = Lv(Lv~=0);
Z = zeros(n,1);
r = 0; % Initialize iterations counter
tol = 1e-6; % Threshold Tolerance
while  any(B > tol)
    % Initializition and Calculation for Q
    for i = 1:m
        Q(i,1) = A(i,1)/norm(A(:,1));
    end
    % Apply Gram-Schmidt algorithm
    for i = 2:n
        S = zeros(m,1);
        for k = 1:i-1
            S = S + (A(:,i)'*Q(:,k))*Q(:,k);
        end
        Z = A(:,i) - S;
        for j = 1:m
            Q(j,i) = Z(j)/norm(Z);
        end
    end
    % Calculation for R
    for i = 1:m
        for j = i:n
            R(i,j) = Q(:,i)'*A(:,j);
        end
    end
    A = R*Q;% Set A=RQ
    % Check if threshold is achieved for all values
    L = tril(A,-1);Lv = L(:);
    B = abs(Lv(Lv~=0));
    % Increament the iteration
    r =r+1; 
end
fprintf('Eigenvalues using QR method without shifts\n')
Eigenvalues = diag(A)
fprintf('Iterations = %g \n',r)
fprintf('\rEigenvalues obtained using built-in matlab function eig()\n')
Matlab_Eigenvalues = eig(A)
