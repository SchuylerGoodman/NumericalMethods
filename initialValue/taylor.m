function [x fx fdx] = taylor(fty, lowerLimit, upperLimit, initialCondition, intervals)
%TAYLOR Estimate value of initial value problem at the given upper limit.
%
%       Inputs:
%               fty - Function handle or cell array of function handle 
%                   representations of differential equation being estimated
%                   and its derivatives.  Order determined by length(fty)
%               lowerLimit - Lower bound for t
%               upperLimit - Upper bound for t, value at which to estimate
%                   y(t)
%               initialCondition - Initial y(lowerLimit) = initialCondition
%               intervals - Number of desired mesh points - 1
%
%       Outputs:
%               x - Row vector of mesh points used
%               fx - Row vector of estimated values of y(t) at mesh points
%               fdx - Row vector of values of y'(t) using estimated values
%                   of y(t) stored in fx.
%
%       NOTE:
%               x, fx, and fdx can be easily used with HERMITE, CUBIC_NAT,
%               or CUBIC_CLAMP.
%

len = length(fty);
assert(len ~= 0, 'Must input a function to estimate.');
if len == 1
    if isa(fty, 'function_handle')
        fty = {fty};
    end
end

%if len < order
%    syms t y;
%    fty = [fty cell(1, order - len)];
%   for d = len + 1 : order
%        fs = sym(fty{len});
%        fsd = diff(fs, t);
%        fty{d} = matlabFunction(fsd, 'vars', [t y]);
%    end
%end

x = zeros(1, intervals + 1);
fx = x;
fdx = x;

step = ( upperLimit - lowerLimit ) / intervals;
current = lowerLimit;
value = initialCondition;
x(1) = current;
fx(1) = value;
fdx(1) = feval(fty{1}, x(1), fx(1));
disp('i     t       estimate');

for i = 2 : intervals + 1
    value = value + step * tayloreval(fty, current, value, step);
    current = current + step;
    x(i) = current;
    fx(i) = value;
    fdx(i) = feval(fty{1}, x(i), fx(i));
    fprintf('%d     %3.2f    %6.5f\n', i - 1, current, value);
end

end

function value = tayloreval(fty, t, y, step)

value = 0;
for i = 1 : length(fty)
    assert(isa(fty{1}, 'function_handle'), 'Differential equations must be of type funtion_handle.');
    currentStep = step^( i - 1 ) / factorial(i);
    value = value + currentStep * feval(fty{i}, t, y);
end

end

