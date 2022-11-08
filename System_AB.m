function [A,b] = System_AB(n,System)
% System_AB(n,Matrix) generates a system of n linear equations
%
% INPUTS: n = number of linear equations
%         System = 'a' or 'b' selects the system type
%
% OUTPUTS: A = The System Matrix of the selected system
%          b = The output vector of the selected system
%
A = zeros(n,n); % Initialize the System Matrix with zeros
b = zeros(n,1); % Initialize the System output matrix with zeros

if  System~='a'&&System~='b'% Set Default to case 'a'
    fprintf('System type selected not available.System of case (a) is used by Default\n')
    System='a';
end

%% System Case 'a'
if System=='a'
    % Apply the Given formulas
    for i=1:n
        for j=1:n
            if i==j
                A(i,j)= 7;
            elseif i==j-1 || i==j+1
                A(i,j) = -2;
            end
        end
        % Generate the b vector
        b(i) = -3+0.5*i;
    end
end

%% System Case 'b'
if System=='b'
    % Apply the given formulas
    for i=1:n
        for j=1:n
            A(i,j)=3/(7*(i+j+1));
        end
        if mod(i,2) == 1 % for i is odd
            b(i) = 9/(7*i);
        else
            b(i) = 0; % for i is even
        end
    end
end
end