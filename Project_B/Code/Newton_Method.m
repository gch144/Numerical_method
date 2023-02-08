%% Netwon Function
function [X,it] = Newton_Method(f, df, x, tol)
maxCount = 50; % Limit iterations to a maximum of 100
counter = 1 ; % Initalize the iterations counter
while abs(f(x))>tol ||maxCount<maxCount
    if df(x)==0
        disp('Division by zero encountered')
        break
    end
    x = x - f(x)/ df(x);    % Newton Raphson Formula
    X(counter) = x;         % Store approximate roots in all iterations
    counter = counter+1;    % Increament the iteration counter
end
it = counter-1; % Total Number of Iterations
end