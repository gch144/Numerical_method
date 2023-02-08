%% Muller's Method 1
function X = MM1(P, x0, x1, x2,tol)
% Evaluate the Intial points
fx0 = polyval(P,x0);
fx1 = polyval(P,x1);
fx2 = polyval(P,x2);

iter = 1; % Start Iterations at 1
I_max = 100; % Maximum number of iterations allowed 
while abs(fx2) > tol && iter <100 I_max;
    % Implement MM1 formulas
    H1 = x1 - x0;   D1 = (fx1 - fx0) / H1;
    H2 = x2 - x1;   D2 = (fx2 - fx1) / H2;
    
    % Calculate interpolation polynomial constants  
    C = (D2 - D1)/(H1 + H2);
    B = D2 + H2*C;
    A = fx2;
    Den = sqrt(B*B + 4*A*C);    
    x0 = x1; x1 = x2;
    
    % Get the minimum of the roots
    if abs(B + Den)> abs(B - Den)
        x2 = x2 - 2*A/(B + Den);
    else
        x2 = x2 - 2*A/(B - Den);
    end
    
    % Update the variables
    X(iter) = x2;
    fx0 = fx1;fx1 = fx2;
    fx2 = polyval(P,x2);
    iter =iter +1;
end
end