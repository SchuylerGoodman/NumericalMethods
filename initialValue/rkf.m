function [t, yt, ydt] = rkf(fty, lowerLimit, upperLimit, initialCondition, tolerance, stepmin, stepmax)
%RKF Estimate value of initial value problem at the given upper limit
%             using Runge-Kutta-Fehlberg Method.
%
%       Inputs:
%               fty - f(t, y), or y'.  The differential equation whose
%                   solution is being estimated.
%               lowerLimit - Lower bound for t
%               upperLimit - Upper bound for t, value at which to estimate
%                   y(t)
%               initialCondition - Initial y(lowerLimit) = initialCondition
%               tolerance - The tolerance to which the results should be
%                   calculated.
%               stepmin - Minimum step size per iteration
%               stepmax - Maximum step size per iteration
%
%       Outputs:
%               t - Row vector of mesh points used
%               yt - Row vector of estimated values of y(t) at mesh points
%               ydt - Row vector of values of y'(t) using estimated values
%                   of y(t) stored in fx.
%
%       NOTE:
%               t, yt, and ydt can be easily used with HERMITE, CUBIC_NAT, and
%                   CUBIC_CLAMP to interpolate other values of y(t).
%


step = stepmax;
t = lowerLimit;
yt = initialCondition;
ydt = feval(fty, lowerLimit, initialCondition);
disp('i     t          estimate');

flag = 1;
i = 1;
while flag
    k1 = step * feval(fty, t(i), yt(i));
    k2 = step * feval(fty, t(i) + 1 / 4 * step, yt(i) + 1 / 4 * k1);
    k3 = step * feval(fty, t(i) + 3 / 8 * step, yt(i) + 3 / 32 * k1 + 9 / 32 * k2);
    k4 = step * feval(fty, t(i) + 12 / 13 * step, yt(i) + 1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3);
    k5 = step * feval(fty, t(i) + step, yt(i) + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4);
    k6 = step * feval(fty, t(i) + 1 / 2 * step, yt(i) - 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5);
    R = 1 / step * abs(1 / 360 * k1 - 128 / 4275 * k3 - 2197 / 75240 * k4 + 1 / 50 * k5 + 2 / 55 * k6); % R = 1/h|w~i+1-wi+1|
    
    if R <= tolerance
        w = yt(i) + 25 / 216 * k1 + 1408 / 2565 * k3 + 2197 / 4104 * k4 - 1 / 5 * k5;
        tnext = t(i) + step;
        t = [t tnext];
        yt = [yt w];
        fprintf('%d     %6.5f    %6.5f\n', i - 1, t(i), yt(i));
        i = i + 1;
        ydt = feval(fty, t(i), yt(i));
    end
    
    del = 0.84 * (tolerance / R)^(1/4);
    if del <= 0.1
        step = 0.1 * step;
    elseif del >= 4
        step = 4 * step;
    else
        step = del * step;
    end
    
    if step > stepmax
        step = stepmax;
    end
    
    if t(i) >= upperLimit
        flag = 0;
    elseif t(i) + step > upperLimit
        step = upperLimit - t(i);
    elseif step < stepmin
        flag = 0;
        disp('Minimum h exceeded');
    end
        
end
fprintf('%d     %6.5f    %6.5f\n', i - 1, t(i), yt(i));

end
