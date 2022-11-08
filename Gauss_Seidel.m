function [x, Rnorms,iter_K] = Gauss_Seidel(A, b, x)
% Gauss_Seidel(A,b, x) solves a system of linear equations by applying
% Gauss Seidel Iterative method
%
% INPUTS A = System matrix/Coefficient matrix
%        b = System Output column vector
%        x = Initial solution vector (same size as vector b)
%
% OUTPUTS: x = final solution of the system
%     Rnorms = a vector of error norms
%     iter_K = a list of iterations used to achieve the solution
%
n = size(A,1);
xnew = x;
Rnorms = [];% Set a vector of Error norms with initial zeros
iter_K = [];% Variable to store the iterations
iter = 1;% Iteration counter

% Applying the Gauss_Seidel iterative algorithm
while norm(A*x - b)>=1e-10% tolerance set to error less than 1e-10
    for i = 1:n
        x_temp = 0;
        for j = 1:n
            if j ~= i
                % Update solution
                x_temp = x_temp + xnew(j) * A(i,j);
            end
        end
        %Calculate a new solution
        xnew(i) = 1/A(i,i)*(b(i) - x_temp);
    end
    x = xnew;               % save the new solution 
    % Calculate for Error norm
    Rnorms(iter) = norm(A*x - b);
    % Store the current iteration value in the variable iter_K
    iter_K(iter) = iter; 
    % Increament the Iteration Counter by 1
    iter = iter+1; 
    % Limit Iteration to a maximum of 10000
    if iter>10000,break,end 
end
end