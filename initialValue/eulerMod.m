function [t, yt, ydt] = eulerMod(fty, lowerLimit, upperLimit, initialCondition, intervals)
%EULERMOD Estimate value of initial value problem at the given upper limit
%             using Modified Euler Method.
%
%       Inputs:
%               fty - f(t, y), or y'.  The differential equation whose
%                   solution is being estimated.
%               lowerLimit - Lower bound for t
%               upperLimit - Upper bound for t, value at which to estimate
%                   y(t)
%               initialCondition - Initial y(lowerLimit) = initialCondition
%               intervals - Number of desired mesh points - 1
%
%       Outputs:
%               t - Row vector of mesh points used
%               yt - Row vector of estimated values of y(t) at mesh points
%               ydt - Row vector of values of y'(t) using estimated values
%                   of y(t) stored in fx.
%
%       NOTE:
%               t and yt can be easily used with HERMITE, CUBIC_NAT, and
%                   CUBIC_CLAMP to interpolate other values of y(t).
%

t = zeros(1, intervals + 1);
yt = t;
ydt = t;

step = ( upperLimit - lowerLimit ) / intervals;
t(1) = lowerLimit;
yt(1) = initialCondition;
ydt(1) = feval(fty, lowerLimit, initialCondition);
disp('i     t       estimate');

for i = 2 : intervals + 1
    t(i) = t(i - 1) + step;
    % yi = yi-1 + h / 2 * (f(ti-1, yi-1) + f(ti, yi-1 + h * f(ti-1, yi-1)))
    a = feval(fty, t(i - 1), yt(i - 1));
    b = yt(i - 1) + step * a;
    c = a + feval(fty, t(i), b);
    yt(i) = yt(i - 1) + step / 2 * c;
    ydt(i) = feval(fty, t(i), yt(i));
    fprintf('%d     %3.2f    %6.5f\n', i - 1, t(i), yt(i));
end

end
