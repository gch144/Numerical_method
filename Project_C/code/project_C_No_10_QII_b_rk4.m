%% Numerical Methods, PROJECT C No. 10
% Problem II: Determine the trajectory of the motion on the interval [0,
% 20] for the following initial conditions: x1(0) = 0.002, x2(0) = 0.02.
% Evaluate the solution using Runge-Kutta method of 4th order (RK4) with a
% variable step size automatically adjusted by the algorithm, making error
% estimation according to the step-doubling rule.
clc;close all;clear all
interval = [0 20];      % Time interval
x0 =[0.002,0.02];       % x0 = [x1(0), x2(0)] 
 
% Define the Equations
f = @(t,x) [x(2)+x(1)*(0.5-x(1)^2-x(2)^2);
    -x(1)+x(2)*(0.5-x(1)^2-x(2)^2)];
 
% Record execution time
tStart = tic;
[X, T,h, err, iter] = RK4_variable_h(f,interval, x0);
tEnd = toc(tStart);
display(iter,'Number of iterations')
display(tEnd ,'Time of execution')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ode45 solver
[t45,x45] = ode45(f,interval,x0);
 
%% Plot Problem solution versus time--
figure()% Create a new figure window
plot(X(:,1),X(:,2),'-b','linewidth',1.5)
title('Solution curves x2 versus x1 using variable h RK4')
xlabel('x1');ylabel('x2') % Add axis labels
%% Plot problem solution vs time 
figure()    % Create a new figure window
plot(T,X,'linewidth',1.5)   % Plot problem solution vs time 
title('Problem solution versus time using variable h RK4')
xlabel('time');ylabel('x1,x2')      % Add axis labels
legend('x1','x2','location','best') % Add legend
 
%% Plot step size vs time 
figure()% Create a new figure window
plot(T,h,'linewidth',1.5)% Plot step size vs time 
title('Step size vs. time')% Add title
xlabel('time');ylabel('Step size,h') % Add axis labels
 
%% Plot error estimate vs time 
figure()% Create a new figure window
plot(T(1:end-1),err,'linewidth',1.5)% Plot step size vs time 
title('Error estimate vs. time') % Add title
xlabel('time');ylabel('Error')  % Add axis labels
legend('error x1','error x2')   % Add legend
 
%% Compare adaptive RK4 and ode45
figure()% Create a new figure window
hold on
plot(x45(:,1),x45(:,2),'-r','linewidth',1.5)
plot(X(:,1),X(:,2),'--b','linewidth',1.5)
title('Comparison of solution of adaptive RK4 and ode45')
xlabel('x1');ylabel('x2') % Add axis labels
legend('Adaptive RK4','ode45') % Add legend
hold off

%% Compare Adaptive RK4 and ode45
figure()% Create a new figure window
hold on
plot(t45, x45,'-r','linewidth',1.5)
plot(T,X,'--b','linewidth',1.5)
title('Comparison of Problem solution versus time')
xlabel('time');ylabel('x1,x2') % Add axis labels
legend('Adaptive RK4 x1','Adaptive RK4 x2','ode45 x1','ode45 x2',...
    'location','best') % Add legends
hold off
