% This scripts implements a program to find macheps in the MATLAB 
% environment on a lab computer or personal computer.
% Clear and clean the Matlab workspace
clc; clear all; close all
macheps = 1.0;
% The condition to stop is: 1+macheps <= 1, 
% thus the condition to continue is 1+macheps > 1
while (1+macheps) > 1
    % Decrease macheps by a factor of 2
    macheps = macheps/2;
end
fprintf('Calculated Macheps value\n macheps = %g\n',macheps*2)
fprintf('\rValue of Matlab''s in-built epsilon value\n')
fprintf('Matlab''s macheps= %g\n',eps)