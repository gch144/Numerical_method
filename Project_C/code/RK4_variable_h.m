function [X,t, h, X_err, counter] = RK4_variable_h(f,interval, y0)
beta = 0.9;             % Safety factor
abs_err = 1e-6;         % Absolute error tolerance
rel_err = 1e-6;         % Relative error tolerance
hmin = 1e-6;            % Set the Minimum step size
cf = 4;                 % Converge order of RK4 method
dh = 0.001;             % Set the initial step size
X(1, :) = y0;           % Record the initial solution
X_err(1,:) = [0, 0];    % Record starting error as zero
t = interval(1);        % Set starting time of the solution
j = 1;                  % usinh symbol  as steps counter
counter = 0;            % Starting counting iteration from 1 (default=0)
h(1) = dh;              % Record first step size

% Apply step-doubling rule to vary the step sizes
while t(j)~= interval(2)
    [ ~, YY ] = RK4_h_constant( f, t(j), t(j)+dh, X(j, :), dh);%Step-size = h
    x1 = YY(end,:);
    dhh = dh/2; %Step-size = h/2
    [ ~, YY ] = RK4_h_constant( f, t(j), t(j)+dhh, X(j, :), dhh);
    x2 = YY(end,:);
    
    delta = abs((x2 - x1))./(2.^cf - 1);
    Error_estimate = abs(x2)*rel_err + abs_err; % Error estimate
    X_err(j, :) = delta; % Record Error
    da = (min(Error_estimate./delta)).^(1/(cf+1)); % Correct step size
    
    B = beta*da;% Calculate Safety factor
    htest = B*dh; % Calculate new step size
    if beta*da>= 1
        t(j+1) = t(j) + dh/2;   % Update time
        X(j+1, :) = x2;         % Update solution
        dh = min([htest, 2*dh, interval(2)- t(j+1)]); % Step size to be used
        j = j + 1; % Increament step counter
        h(j) = dh; % Record the step size
    else
        if htest < hmin
            disp('Error!:Minimum stepp size reached!.');
            break;
        else
            dh = htest;
        end
    end
    counter = counter + 1;% Update iteration counter
end
end
