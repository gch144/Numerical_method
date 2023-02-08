%% Numerical Methods, PROJECT C No. 10
% Problem II:
% Determine the trajectory of the motion on the interval [0, 20] for the
% following initial conditions: x1(0) = 0.002, x2(0) = 0.02. Evaluate the
% solution using: a) Runge-Kutta method of 4th order (RK4)
clc;close all;clear all
% Define the Equations

f = @(t,x) [x(2)+x(1)*(0.5-x(1)^2-x(2)^2);
    -x(1)+x(2)*(0.5-x(1)^2-x(2)^2)];
% Define the initial conditions
t0 = 0;             % Start Time
tf = 20;            % End time
interval = [t0,tf]; % Time interval
x0 =[0.002,0.02];   % x0 = [x1(0), x2(0)]

%% Solution
% Evaluate the solution using different constant step-sizes
ode45_sol = ode45(f,interval,x0);
dh = 0.05:0.1:1.1;
sol_error = zeros(length(dh),2);
for i=1:length(dh)
    [ rk4_t, rk4_sol ] = RK4_h_constant(f, t0, tf, x0, dh(i));
    ode45_sol_x = deval(ode45_sol,rk4_t);
    abs_error = abs(rk4_sol-ode45_sol_x');
    sol_error(i,1)=max(abs_error(:,1));
    sol_error(i,2)=max(abs_error(:,2));
end
% Compare errors for each step size.
table(dh(:),sol_error(:,1),sol_error(:,2),...
    'variablenames',{'h','del_x1','del_x2'})
% Plot the error trend with increasing step size, h
figure()                                % Create a new figure window
plot(dh,sol_error,'linewidth',1.5)      % Add the plot on the figure
grid on                                 % Add grid to the plot
xlabel('Step size, h');ylabel('Error')  % Add axis labels
legend('x1 error','x2 error','location','best')% Add a legend
title('Solution errors of RK4 compared to solution of ode45')  % Add title

%% Plot Solution curves x2 versus x1
[t45,ode45_sol_x] = ode45(f,interval,x0);
[ Trk4_optimal, Xrk4_optimal ] = RK4_h_constant(f, t0, tf, x0, 0.1);
[ Trk4_large, Xrk4_large ] = RK4_h_constant(f, t0, tf, x0, 0.75);
figure() % Create a new figure window
hold on  % holds the current plot and all axis properties
plot(Xrk4_optimal(:,1),Xrk4_optimal(:,2),'-b','linewidth',1.5)
plot(Xrk4_large(:,1),Xrk4_large(:,2),'--r','linewidth',1.5)
title('Solution curves x2 versus x1 using RK4') % Add title
xlabel('x1');ylabel('x2')                       % Add axis labels
legend('optimal h','large h','location','best') % Add legend
grid on                                         % Add grid lines
hold off
%% Plot Problem solution versus time
figure() % Create a new figure window
hold on  % holds the current plot and all axis properties
plot(Trk4_optimal,Xrk4_optimal,'-r','linewidth',1.5)
plot(Trk4_large, Xrk4_large,'--k','linewidth',1.5)
plot(t45, ode45_sol_x,'--g','linewidth',1.5)
title('Problem solution versus time using RK4') % Add title
xlabel('time');ylabel('x1,x2')                  % Add axis labels
legend('x1(optimal h)','x2(optimal h)','x1(large h)','x2(large h)',...
    'x1(ode45)','x2(ode5)','location','best')   % Add legend
grid on                                         % Add grid lines
hold off
%% Compare the results with the ones obtained using ode45.
figure() % Create a new figure window
hold on  % holds the current plot and all axis properties
plot(Xrk4_optimal(:,1),Xrk4_optimal(:,2),'-r','linewidth',1.5)
plot(ode45_sol_x(:,1),ode45_sol_x(:,2),'--g','linewidth',1.5)
title('Comparison of solution of RK4 and ode45')    % Add title
xlabel('x1');ylabel('x2')                           % Add axis labels
legend('RK4','ode45','location','best')             % Add legend
grid on                                             % Add grid lines
hold off