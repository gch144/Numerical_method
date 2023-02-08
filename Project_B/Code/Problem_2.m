% Numerical Methods, project B No. 10
% Find all (real and complex) roots of the polynomial
% f(x) = a4x4+ a3x3+ a2x2+a1x+a0, [a4 a3 a2 a1 a0] = [1 7 -5 8 9]
% a) the Muller method 1,
% b) the Muller method 2,
% c) the Newton's method.
%
% Clear the work space
clc;clear;close all;
format short;

% Define the function and its derivatives.
P = [1 7 -5 8 9];           % Coefficients of the polynomial
dP = polyder(P);            % Coefficients of the 1st derivative
dP2 = polyder(dP);          % Coefficients of the 2nd derivative
f = @(x) polyval(P,x);      % Anonymous function for the polynomial
df = @(x) polyval(dP,x);    % Anonymous function for the 1st derivative
df2 = @(x) polyval(dP2,x);  % Anonymous function for the 2nd derivative

% Define other constants to use in the program
tol = 1e-10; % Set the toleance
mm1_roots = zeros(4,1);
mm2_roots = zeros(4,1);
newt_roots = zeros(2,1);

% Visualize the function revealing all x-intercepts 
figure() % Figure 1
fplot(f,[-8.5 ,3.5])
yline(0);xline(0);grid on;
xlabel('x');ylabel('f(x)');grid on 
title('f(x)= x^4 + 7x^3 - 5x^2 + 8x + 9 ')


%% Solve for all Roots using MM1
% Finding first root
IP1 = [-6,-7,-8];% First set of INITIAL POINTS
mm1_X1 = MM1(P, IP1(1), IP1(2), IP1(3), tol);
mm1_roots(1) = mm1_X1(end); % Store the first Root

% Finding second root
IP2 = [-2,1,2];% Second set of INITIAL POINTS
mm1_X2 = MM1(P, IP2(1), IP2(2), IP2(3), tol);
mm1_roots(2) = mm1_X2(end); % Store the second Root

% Finding third root
IP3 = [1+0.5i,1+1i,1+1.5i];% Third set of INITIAL POINTS
mm1_X3 = MM1(P, IP3(1), IP3(2), IP3(3), tol);
mm1_roots(3) = mm1_X3(end); % Store the third Root

% Finding fourth root
IP4 = [1-0.5i,1-1i,1-1.5i];% Fourth set of INITIAL POINTS
mm1_X4 = MM1(P, IP4(1), IP4(2), IP4(3), tol);
mm1_roots(4) = mm1_X4(end); % Store the fourth Root

% Display Results 
% display(mm1_roots, 'Roots by MM1')

% Plot the results
figure;hold on;
fplot(f,[-8.5 ,3.5]);
realRoots = mm1_roots(abs(imag(mm1_roots)) <1e-10);
rr = real(realRoots);frr = f(rr);
plot(IP1,f(IP1),'sr','MarkerFaceColor','r')
plot(IP2,f(IP2),'sg','MarkerFaceColor','g')
plot(rr,frr,'*r')
yline(0);xline(0)
grid on;xlabel('x');ylabel('f(x)');
title('Real roots with MM1');
legend('f(x)','1st Intial points','2nd Intial points','roots',...
    'location','southeast');

%% Solve for all Roots using MM2)
IG = [-6, -2,1+0.5i,1-0.5i];% mm2 Initial Points

% Finding first root
mm2_X1 = MM2(P, dP, dP2, IG(1), tol);
mm2_roots(1) = mm2_X1(end); % Store the first Root

% Finding second root
mm2_X2 = MM2(P, dP, dP2, IG(2), tol);
mm2_roots(2) = mm2_X2(end); % Store the second Root

% Finding third root
mm2_X3 = MM2(P, dP, dP2, IG(3), tol);
mm2_roots(3) = mm2_X3(end); % Store the third Root

% Finding fourth root
mm2_X4 = MM2(P, dP, dP2, IG(4), tol);
mm2_roots(4) = mm2_X4(end); % Store the fourth Root

% Plot the results
figure;hold on; 
fplot(f,[-8.5 ,3.5]);
realRoots2 = mm2_roots(abs(imag(mm2_roots)) <1e-10);
rr2 = real(realRoots2);frr2 = f(rr2);
plot(IG(1),f(IG(1)),'sr','MarkerFaceColor','r')
plot(IG(2),f(IG(2)),'sg','MarkerFaceColor','g')
plot(rr2,frr2,'*r')
yline(0);xline(0)
grid on;xlabel('x');ylabel('f(x)');
title('Real roots with MM2');
legend('f(x)','1st Intial point','2nd Intial point','roots',...
    'location','southeast');

%% APPLY NEWTON METHOD
% Finding first root
newt_X1 = Newton_Method(f, df, IG(1), tol);
newt_roots(1) = newt_X1(end); % Store the first Root

% Finding second root
newt_X2 = Newton_Method(f, df, IG(2), tol);
newt_roots(2) = newt_X2(end); % Store the first Root

% Plot the results
figure;hold on;
fplot(f,[-8.5 ,3.5]);
plot(IG(1),f(IG(1)),'sr','MarkerFaceColor','r')
plot(IG(2),f(IG(2)),'sg','MarkerFaceColor','g')
plot(newt_roots,f(newt_roots),'*r')
yline(0);xline(0)
grid on;xlabel('x');ylabel('f(x)');
title('Real roots with NEWTON METHOD');
legend('f(x)','1st Intial point','2nd Intial point','roots',...
    'location','southeast');

%% Display and Compare the methods
fprintf('Roots by MM1\n')
nMM1 = [length(mm1_X1),length(mm1_X2),length(mm1_X3),length(mm1_X4)];
MM1sol = [mm1_roots(:),f(mm1_roots(:)),nMM1(:)];
fprintf('r1 = %7.4f%+7.4fi      f(r1) = %7.5e    Iter = %g\n',...
    real(MM1sol(1,1)),imag(MM1sol(1,1)),MM1sol(1,2),MM1sol(1,3))
fprintf('r2 = %7.4f%+7.4fi      f(r2) = %7.5e    Iter = %g\n',...
    real(MM1sol(2,1)),imag(MM1sol(2,1)),MM1sol(2,2),MM1sol(2,3))
fprintf('r3 = %7.4f%+7.4fi      f(r3) = %7.5e    Iter = %g\n',...
    real(MM1sol(3,1)),imag(MM1sol(3,1)),MM1sol(3,2),MM1sol(3,3))
fprintf('r4 = %7.4f%+7.4fi      f(r4) = %7.5e    Iter = %g\n',...
    real(MM1sol(4,1)),imag(MM1sol(4,1)),MM1sol(4,2),MM1sol(4,3))

fprintf('\r\rRoots by MM2\n')
nMM2 = [length(mm2_X1),length(mm2_X2),length(mm2_X3),length(mm2_X4)];
MM2sol = [mm2_roots(:),f(mm2_roots(:)),nMM2(:)];
fprintf('r1 = %7.4f%+7.4fi      f(r1) = %7.5e    Iter = %g\n',...
    real(MM2sol(1,1)),imag(MM2sol(1,1)),MM2sol(1,2),MM2sol(1,3))
fprintf('r2 = %7.4f%+7.4fi      f(r2) = %7.5e    Iter = %g\n',...
    real(MM2sol(2,1)),imag(MM2sol(2,1)),MM2sol(2,2),MM2sol(2,3))
fprintf('r3 = %7.4f%+7.4fi      f(r3) = %7.5e    Iter = %g\n',...
    real(MM2sol(3,1)),imag(MM2sol(3,1)),MM2sol(3,2),MM2sol(3,3))
fprintf('r4 = %7.4f%+7.4fi      f(r4) = %7.5e    Iter = %g\n',...
    real(MM2sol(4,1)),imag(MM2sol(4,1)),MM2sol(4,2),MM2sol(4,3))

fprintf('\r\rRoots by Newton Method\n')
nNewt = [length(newt_X1),length(newt_X2)];
newtsol = [newt_roots(:),f(newt_roots(:)),nNewt(:)];
fprintf('r1 = %7.4f     f(r1) = %7.5e    Iter = %g\n',...
    newtsol(1,1),newtsol(1,2),newtsol(1,3))
fprintf('r2 = %7.4f     f(r2) = %7.5e    Iter = %g\n',...
    newtsol(2,1),newtsol(2,2),newtsol(2,3))

% A comparison of the results containing a table with all successive
% iteration points for the first root
%************************************************************************
fprintf('\r\r')
iter_mm1 = [1:length(mm1_X1)]';
roots_mm1 = mm1_X1(:);
funValue_mm1 = f(mm1_X1(:));
MM1_Table = table(iter_mm1,roots_mm1,funValue_mm1); 
disp(MM1_Table)
%************************************************************************
iter_mm2 = [1:length(mm2_X1)]';
roots_mm2 = mm2_X1(:);
funValue_mm2 = f(mm2_X1(:));
MM2_Table = table(iter_mm2,roots_mm2,funValue_mm2); 
disp(MM2_Table)

%************************************************************************
iter_newt = [1:length(newt_X1)]';
roots_newt = newt_X1(:);
funValue_newt = f(newt_X1(:));
newt_Table = table(iter_newt,roots_newt,funValue_newt); 
disp(newt_Table)