% Numerical Methods, project B No. 10
% Find all (real and complex) roots of the polynomial
% f(x) = a4x4+ a3x3+ a2x2+a1x+a0, [a4 a3 a2 a1 a0] = [1 7 -5 8 9]
% a) Laguerre’s method ,
% b) the Muller method 2,
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
L_roots = zeros(4,1);
mm2_roots = zeros(4,1);
IG = [-6, -2,1+0.5i,1-0.5i];% Initial guess Points

% Visualize the function revealing all x-intercepts 
figure() % Figure 1
fplot(f,[-8.5 ,3.5])
yline(0);xline(0);grid on;
xlabel('x');ylabel('f(x)');grid on 
title('f(x)= x^4 + 7x^3 - 5x^2 + 8x + 9 ')


%% Solve for all Roots using Laguerre Method)
% Finding first root
L_X1 = laguerre_Method(P, dP, dP2, IG(1), tol);
L_roots(1) = L_X1(end); % Store the first Root

% Finding Second root
L_X2 = laguerre_Method(P, dP, dP2, IG(2), tol);
L_roots(2) = L_X2(end); % Store the second Root

% Finding Third root
L_X3 = laguerre_Method(P, dP, dP2, IG(3), tol);
L_roots(3) = L_X3(end); % Store the Third Root

% Finding Third root
L_X4 = laguerre_Method(P, dP, dP2, IG(4), tol);
L_roots(4) = L_X4(end); % Store the Third Root

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Plot the results
% Plot the results
figure;hold on; 
fplot(f,[-8.5 ,3.5]);
realRoots_L = L_roots(abs(imag(L_roots)) <1e-10);
Lr = real(realRoots_L);fLr = f(Lr);
plot(IG(1),f(IG(1)),'sr','MarkerFaceColor','r')
plot(IG(2),f(IG(2)),'sg','MarkerFaceColor','g')
plot(Lr,fLr,'*r')
yline(0);xline(0)
grid on;xlabel('x');ylabel('f(x)');
title('Real root with Laguerre’s method');
legend('f(x)','1st Intial point','2nd Intial point','roots',...
    'location','southeast');

%% Solve for all Roots using MM2)
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
realRoots_L = mm2_roots(abs(imag(mm2_roots)) <1e-10);
rr2 = real(realRoots_L);frr2 = f(rr2);
plot(IG(1),f(IG(1)),'sr','MarkerFaceColor','r')
plot(IG(2),f(IG(2)),'sg','MarkerFaceColor','g')
plot(rr2,frr2,'*r')
yline(0);xline(0)
grid on;xlabel('x');ylabel('f(x)');
title('Real root with MM2');
legend('f(x)','1st Intial point','2nd Intial point','roots',...
    'location','southeast');

%% Display and Compare the methods
fprintf('Roots by Laguerre’s method\n')
nL = [length(L_X1),length(L_X2),length(L_X3),length(L_X4)];
L_sol = [L_roots(:),f(L_roots(:)),nL(:)];
fprintf('r1 = %7.4f%+7.4fi      f(r1) = %7.5e    Iter = %g\n',...
    real(L_sol(1,1)),imag(L_sol(1,1)),L_sol(1,2),L_sol(1,3))
fprintf('r2 = %7.4f%+7.4fi      f(r2) = %7.5e    Iter = %g\n',...
    real(L_sol(2,1)),imag(L_sol(2,1)),L_sol(2,2),L_sol(2,3))
fprintf('r3 = %7.4f%+7.4fi      f(r3) = %7.5e    Iter = %g\n',...
    real(L_sol(3,1)),imag(L_sol(3,1)),L_sol(3,2),L_sol(3,3))
fprintf('r4 = %7.4f%+7.4fi      f(r4) = %7.5e    Iter = %g\n',...
    real(L_sol(4,1)),imag(L_sol(4,1)),L_sol(4,2),L_sol(4,3))

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

% A comparison of the results containing a table with all successive
% iteration points for the first root
%************************************************************************
fprintf('\r\r')
iter_Laguerre = [1:length(L_X1)]';
roots_Laguerre = L_X1(:);
funValue_Laguerre = f(L_X1(:));
Laguerre_Table = table(iter_Laguerre,roots_Laguerre,funValue_Laguerre); 
disp(Laguerre_Table)
%************************************************************************
iter_mm2 = [1:length(mm2_X1)]';
roots_mm2 = mm2_X1(:);
funValue_mm2 = f(mm2_X1(:));
MM2_Table = table(iter_mm2,roots_mm2,funValue_mm2); 
disp(MM2_Table)

