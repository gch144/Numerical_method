function [ T, Y ] = RK4_h_constant( f, t0, tf, y0, h)
T = (t0 : h : tf)'; % Time vector
Y = zeros(length(T),length(y0)); % Initialize solution vector
Y(1,1) = y0(1); % Append initial conditions to the solution vector
Y(1,2) = y0(2);% Append initial conditions to the solution vector
% Implementing the Runge-Kutta method of 4th order (RK4) algorithm
for i = 2 : length(T)
    k1 = f(T(i-1),Y(i-1,:))';
    k2 = f(T(i-1)+0.5*h,Y(i-1,:)+0.5*h*k1)';
    k3 = f(T(i-1)+0.5*h,Y(i-1,:)+0.5*h*k2)';
    k4 = f(T(i-1)+h,Y(i-1,:)+h*k3)';
    Y(i,:) = Y(i-1,:) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end
end
