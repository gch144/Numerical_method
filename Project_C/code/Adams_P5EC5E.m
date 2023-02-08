function [T,Y] = Adams_P5EC5E(dyfun, tspan, y0, h)
k = 5; % Use 5 steps in the preictor and collector
% Runge-Kutta method is used to Initialize the first 5 points
ti = tspan(1):h:tspan(2);
[~,yy] = RK4_h_constant(dyfun,ti(1),ti(5),y0,h);
Y(1:5, :)= yy;
T(1:5) = ti(1:5);
 
% Predictor function
yn0 = @(y, f, h) (y + h/720*(1901*f(end,:)-2774*f(end-1,:) + ...
    2616*f(end-2,:) - 1274*f(end-3,:) + 251*f(end-4,:)));
 
% Corrector function
yn = @(y, f, h, fn0) (y + h/1440*(1427*f(end,:) - 798*f(end-1,:) + ...
    482*f(end-2,:) - 173*f(end-3,:) + 27*f(end-4,:)) + h*475/1440*fn0);
 
f = zeros(length(Y), 2); % Initialize f same size as X with zeros
N = (tspan(2) - tspan(1)) / h; % Total number of steps
% Adams PC (P5EC5E)
for i = k+1:N
    T(i) = T(i-1) + h;
    % P-prediction
    P = yn0(Y(i-1, :), f(1:i-1, :), h);
    
    % E-Evaluation
    f(i,:) = dyfun(T(i),P)';
    % C-Correction
    C = yn(Y(i-1, :), f(1:i-1, :), h, f(i, :));
    
    % E-Evaluation
    f(i,:) = dyfun(T(i),C)';
    
    Y(i, :) = C;
end
T = T(:);
end
