%% Function for Muller's Method 2
function X= MM2(P, dP, dP2, x, tol)
% This method requires a single intial point
C = polyval(P,x);
count = 1;
% Implement muller's method 2
while abs(C) > tol && count<100
    A = 0.5 * polyval(dP2,x);
    B = polyval(dP,x);
    D = sqrt(B * B - 4 * A * C);
    
    % Check and use maximum term of denominator
    if abs(B + D) > abs(B - D)
        x = x - 2 * C / (B + D);
    else
        x = x - 2 * C / (B - D);
    end
    
    X(count) = x;
    C = polyval(P,x);
end
end
