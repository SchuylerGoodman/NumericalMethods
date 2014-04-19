function [x fx] = euler(fty, lowerLimit, upperLimit, initialCondition, intervals)
%EULER Estimate value of initial value problem at the given upper limit.
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
%               x - Row vector of mesh points used
%               fx - Row vector of estimated values of y(t) at mesh points
%
%       NOTE:
%               x and fx can be easily used with CUBIC_NAT and CUBIC_CLAMP
%               to estimate other values of y(t).
%

x = zeros(1, intervals + 1);
fx = x;

step = ( upperLimit - lowerLimit ) / intervals;
current = lowerLimit;
value = initialCondition;
x(1) = current;
fx(1) = value;
disp('i     t       estimate');

for i = 2 : intervals + 1
    value = value + step * feval(fty, current, value);
    current = current + step;
    x(i) = current;
    fx(i) = value;
    fprintf('%d     %3.2f    %6.5f\n', i - 1, current, value);
end

end

