%% Numerical Methods, PROJECT C No. 10
% Problem II:
% Determine the trajectory of the motion on the interval [0, 20] for the
% following initial conditions: x1(0) = 0.002, x2(0) = 0.02. Evaluate the
% solution using: a) Adams PC (P5EC5E)
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
dh = [0.8,0.5,0.4,0.3,0.2,0.1,0.01,0.005,0.001,0.0005];
sol_error = zeros(length(dh),2);
for i=1:length(dh)
    [ T_adams, X_Adams ] = Adams_P5EC5E(f, interval,x0, dh(i));
    ode45_sol_x = deval(ode45_sol,T_adams);
    abs_error = abs(X_Adams-ode45_sol_x');
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
title('Solution errors of Adams PC compared to solution of ode45')% title 
%% Plot Solution curves x2 versus x1
[t45,ode45_sol_x] = ode45(f,interval,x0);
[ Tadams_optimal, Xadams_optimal ] = Adams_P5EC5E(f, interval,x0, 0.001);
[ Tadams_large, Xadams_large ] = Adams_P5EC5E(f, interval,x0, 0.5);
figure() % Create a new figure window
hold on  % holds the current plot and all axis properties
plot(Xadams_optimal(:,1),Xadams_optimal(:,2),'-b','linewidth',1.5)
plot(Xadams_large(:,1),Xadams_large(:,2),'--r','linewidth',1.5)
title('Solution curves x2 versus x1 using Adams PC') % Add title
xlabel('x1');ylabel('x2')                       % Add axis labels
legend('optimal h','large h','location','best') % Add legend
grid on                                         % Add grid lines
hold off
%% Plot Problem solution versus time
figure() % Create a new figure window
hold on  % holds the current plot and all axis properties
plot(Tadams_optimal,Xadams_optimal,'-r','linewidth',1.5)
plot(Tadams_large, Xadams_large,'--k','linewidth',1.5)
plot(t45, ode45_sol_x,'--g','linewidth',1.5)
title('Problem solution versus time using Adams PC') % Add title
xlabel('time');ylabel('x1,x2')                  % Add axis labels
legend('x1(optimal h)','x2(optimal h)','x1(large h)','x2(large h)',...
    'x1(ode45)','x2(ode5)','location','best')   % Add legend
grid on                                         % Add grid lines
hold off
%% Compare the results with the ones obtained using ode45.
figure() % Create a new figure window
hold on  % holds the current plot and all axis properties
plot(Xadams_optimal(:,1),Xadams_optimal(:,2),'-r','linewidth',1.5)
plot(ode45_sol_x(:,1),ode45_sol_x(:,2),'--g','linewidth',1.5)
title('Comparison of solution of Adams PC and ode45')    % Add title
xlabel('x1');ylabel('x2')                           % Add axis labels
legend('Adams PC','ode45','location','best')             % Add legend
grid on                                             % Add grid lines
hold off