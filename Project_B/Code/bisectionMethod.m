%% Bisection function
function [X,it]= bisectionMethod(f, a, b, tol)
N = ((log10(b-a)-log10(tol))/log10(2))-1;
N = ceil(N);
% check validity
if f(a)*f(b)>0
    disp('f(a)*f(b)>0!'),
    disp('No root available in the interval')
    X = NaN;
    return
end
it=0;
X = zeros(1,N);
while it < N
    c = (a+b)/2;    % Callculate new root estimate
    X(it+1) = c;      % Store the new root estimate
    if f(a)*f(c)> 0 % If f(a)*f(b)>0, set a = c
        a = c;
    else            % If f(a)*f(b)<0, set b = c
        b = c;
    end
    it=it+1;        % Increament the iteration
end
end