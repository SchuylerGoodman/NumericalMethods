function [ yt ] = adamsPC( fty, lowerLimit, upperLimit, t, w, steps )
%ADAMSPC Estimate value of initial value problem using the
%           predictor-corrector method with the Adams-Bashfort and
%           Adams-Moulton methods.
%
%       Inputs:
%               fty - f(t, y), or y'.  The differential equation whose
%                   solution is being estimated.
%               lowerLimit - Lower bound for t
%               upperLimit - Upper bound for t, value at which to estimate
%                   y(t)
%               t - Values of t for which y(t) is known. Must be evenly
%                   spaced intervals.
%               w - Values of y(t) corresponding one-to-one to values in t.
%               steps - Integer choice of step.  Possibilities are 3 or 4.
%                   If steps = N, then length(t) == length(w) >= N.
%                   If more values are given in t and w than necessary, the
%                   first N will be used.
%               start - Value of t and w at which to start.
%
%       Outputs:
%               yt - Row vector of exact and approximated values for y(t)
%                   from lowerLimit to upperLimit.
%

assert(nargin > 5, 'Not enough parameters.');
assert(length(t) == length(w), 'Invalid parameters t and w: Must be equal lengths.');
assert(steps > 2 && steps < 5, 'Invalid step size.');
assert(length(t) >= steps, 'Not enough values given for step size.');

step = ( upperLimit - lowerLimit ) / ( length(t) - 1 );
yt = zeros(1, length(w));

j = steps;
yt(1 : j) = w(1 : j);

for i = j : length(t) - 1
    switch steps
        case 3
            wj = yt(i) + step / 12 * ( 23 * feval(fty, t(i), yt(i)) - 16 * feval(fty, t(i - 1), yt(i - 1)) + 5 * feval(fty, t(i - 2), yt(i - 2)));
            wj = yt(i) + step / 12 * ( 5 * feval(fty, t(i) + step, wj) + 8 * feval(fty, t(i), yt(i)) - feval(fty, t(i - 1), yt(i - 1)));
        case 4
            wj = yt(i) + step / 24 * ( 55 * feval(fty, t(i), yt(i)) - 59 * feval(fty, t(i - 1), yt(i - 1)) + 37 * feval(fty, t(i - 2), yt(i - 2)) - 9 * feval(fty, t(i - 3), yt(i - 3)));
            wj = yt(i) + step / 24 * ( 9 * feval(fty, t(i) + step, wj) + 19 * feval(fty, t(i), yt(i)) - 5 * feval(fty, t(i - 1), yt(i - 1)) + feval(fty, t(i - 2), yt(i - 2)));
    end
    yt(i + 1) = wj;
    fprintf('y(%5.3f) = %6.8f\n', t(i + 1), yt(i + 1));
end

end
