function xi = GaussPP(A,b)
% Numerical Methods, project A No. 9
% The indicated method: Gaussian elimination with partial pivoting
% GaussPP(A,b) find the solution of a system described by a system mattrix
% A and an output vector b.
% INPUTS: A = System matrix
%         b = Output column vector
%
% OUTPUTS x = solution of the system
%
% The System matrix must be square
if diff(size(A))~=0, error('Matrix A IS NOT square!'); end

% Form the augmented matrix using A and b
Ab=[A b];

[~,n] = size(A); % Get the Dimensions of matrix A
n1 = n+1;
% Gaussian forward elimination
for k = 1:n-1
    % Partial pivoting
    [~,Pi] = max(abs(Ab(k:n,k))); % Pi is the index of Pivot Row
    Si = Pi+k-1; % Index used to switch rows
    if Si~=k
        % Switch the rows
        Ab([k,Si],:)=Ab([Si,k],:); 
    end
    % Forward elimination
    for Pi = k+1:n 
        % Calculate Multiplication Factor
        MF = Ab(Pi,k)/Ab(k,k); 
        Ab(Pi,k:n1)= Ab(Pi,k:n1)-MF*Ab(k,k:n1);% Elimination
    end
end

% Apply back substitution
xi = zeros(n,1);
% Start by solving the the last value of x
xi(n) = Ab(n,n1)/Ab(n,n); 
for Pi = n-1:-1:1 
    xi(Pi) = (Ab(Pi,n1)-Ab(Pi,Pi+1:n)*xi(Pi+1:n))/Ab(Pi,Pi);
end
end