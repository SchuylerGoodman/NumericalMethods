function [t, yt, ydt] = rkv(fty, lowerLimit, upperLimit, initialCondition, tolerance, stepmin, stepmax)
%RKV Estimate value of initial value problem at the given upper limit
%             using Runge-Kutta-Verner Method.
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
    k2 = step * feval(fty, t(i) + 1 / 6 * step, yt(i) + 1 / 6 * k1);
    k3 = step * feval(fty, t(i) + 4 / 15 * step, yt(i) + 4 / 75 * k1 + 16 / 75 * k2);
    k4 = step * feval(fty, t(i) + 2 / 3 * step, yt(i) + 5 / 6 * k1 - 8 / 3 * k2 + 5 / 2 * k3);
    k5 = step * feval(fty, t(i) + 5 / 6 * step, yt(i) - 165 / 64 * k1 + 55 / 6 * k2 - 425 / 64 * k3 + 85 / 96 * k4);
    k6 = step * feval(fty, t(i) + step, yt(i) + 12 / 5 * k1 - 8 * k2 + 4015 / 612 * k3 - 11 / 36 * k4 + 88 / 255 * k5);
    k7 = step * feval(fty, t(i) + 1 / 15 * step, yt(i) - 8263 / 15000 * k1 + 124 / 75 * k2 - 643 / 680 * k3 - 81 / 250 * k4 + 2484 / 10625 * k5);
    k8 = step * feval(fty, t(i) + step, yt(i) + 3501 / 1720 * k1 - 300 / 43 * k2 + 297275 / 52632 * k3 - 319 / 2322 * k4 + 24068 / 84065 * k5 + 3850 / 26703 * k7);
    R = 1 / step * abs(yt(i) + 3 / 40 * k1 + 875 / 2244 * k3 + 23 / 72 * k4 + 264 / 1955 * k5 + 125 / 11592 * k7 + 43 / 616 * k8 - ( yt(i) + 13 / 160 * k1 + 2375 / 5984 * k3 + 5 / 16 * k4 + 12 / 85 * k5 + 3 / 44 * k6 )); % R = 1/h|w~i+1-wi+1|
    
    if R <= tolerance
        w = yt(i) + 13 / 160 * k1 + 2375 / 5984 * k3 + 5 / 16 * k4 + 12 / 85 * k5 + 3 / 44 * k6;
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
