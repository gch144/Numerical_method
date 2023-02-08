% Numerical Methods, project B No. 10
% Find all zeros of the function
%         f(x0 = 0.4x*cos(5x)-4ln(x+3)
% in the interval [5, 14] using:
% a) the Bisection method,
% b) the Newton's method
%
% Clear the work space
clc;clear;close all;
format long;

% Define the function as anonymous function and plot in the search interval
f = @(x) 0.4.*x.*cos(5.*x)-log(x+3);
figure()
fplot(f,[5,14]);yline(0) % Plot the function in the given interval
grid on;xlabel('x');ylabel('y') % Add axis labels
title('f(x) = 0.4x*cos(5x)-ln(x+3)')

% We notice the function oscillates and therefore, we can use the x-values
% at peaks and troughs to define serch intervals.
x = linspace(5,14,1e4);
[pks, pklocs] = findpeaks(f(x));        % Original Peaks
[troughs,trlocs] = findpeaks(-f(x));    % Original Troughs
pktr = [pklocs(:);trlocs(:)];           % Locations
pktrs = sortrows(pktr);                 % Sorted LocationS
pktrs=pktrs(2:end);                     % Remove first index
s_intervals = x(pktrs)';                % x_values of peaks and troughs
f_intervals = f(s_intervals);           % f_values at the peaks and troughs


%% Roots by the Bisection Method
fprintf('Zeros using the Bisection Method\n')
tol = 1e-10;     % Tolerance
Num_r = length(s_intervals);
Sol_b = zeros(Num_r-1,3);
fprintf('  root      \tf(root) \titeration\n')
for i=1:Num_r-1
    a = s_intervals(i); b = s_intervals(i+1); % Search Interval
    [X,it] = bisectionMethod(f, a, b, tol); % Call Bisection method
    Sol_b(i,1) = X(end); % Root
    Sol_b(i,2) = f(X(end)); % function value at the root
    Sol_b(i,3) = it; % Number of iteration used to obtain the root
    fprintf('%7.4f \t%12.5e\t\t%g\n',Sol_b(i,1),Sol_b(i,2),Sol_b(i,3))
    if i==1
        X_b = X;
    end
end
figure()
plot(x,f(x),'b-');% Plot the function in the given interval
hold on
plot(s_intervals,f_intervals,'ms','MarkerFaceColor','m')
plot(Sol_b(:,1),Sol_b(:,2),'r*'); % Plot the rooots in the given interval
yline(0)
grid on;xlabel('x');ylabel('y') % Add axis labels
title('Roots by Bisection Method')
legend('function','Interval points','zeros','location','northwest')
hold off


%% Roots by the Newton's Method
fprintf('\rZeros using the Newton''s Method\n')
% Determine the delivative of the function
Fsys = sym(f);                 % Convert to symbolic expression
dFsys = diff(Fsys);            % Get the symbolic derivative
df = matlabFunction(dFsys);    % Convert derivative to anonymous function

tol = 1e-10;     % Tolerance
Num_r = length(s_intervals);
Sol_n = zeros(Num_r-1,3);
initial_guess = (s_intervals(1:end-1)+s_intervals(2:end))./2;
fprintf('  root      \tf(root) \titeration\n')
for i=1:Num_r-1
    r = initial_guess(i); % SInitail guess
    [X,it] = Newton_Method(f, df, r, tol); % Call Newtons method
    Sol_n(i,1) = X(end); % Root
    Sol_n(i,2) = f(X(end)); % function value at the root
    Sol_n(i,3) = it; % Number of iteration used to obtain the root
    fprintf('%7.4f \t%12.5e\t\t%g\n',Sol_n(i,1),Sol_n(i,2),Sol_n(i,3))
    if i==1
        X_n = X;
    end
end
figure()
plot(x,f(x),'b-');% Plot the function in the given interval
hold on
plot(initial_guess,f(initial_guess),'ms','MarkerFaceColor','m')
plot(Sol_n(:,1),Sol_n(:,2),'r*'); % Plot the rooots in the given interval
yline(0)
grid on;xlabel('x');ylabel('y') % Add axis labels
title('Roots by Newton''s method')
legend('function','Initial guess','zeros','location','northwest')
hold off

%% Comparison
% A comparison of the results containing a table with all successive
% iteration points for the first root
%**********************************************************************
fprintf('\r\r')
iter_bisection = [1:length(X_b)]';
roots_bisection = X_b(:);
funValue_bisection= f(X_b(:));
B_Table = table(iter_bisection,roots_bisection,funValue_bisection);
disp(B_Table)
%************************************************************************
iter_Newton = [1:length(X_n)]';
roots_Newton = X_n(:);
funValue_Newton = f(X_n(:));
N_Table = table(iter_Newton,roots_Newton,funValue_Newton);
disp(N_Table)