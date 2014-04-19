function [ yt ] = adamsM( fty, lowerLimit, upperLimit, t, w, steps )
%ADAMSM Estimate value of initial value problem using the implicit
%           Adams-Moulton Methods
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
%               steps - Integer choice of step.  Possibilities are 2, 3, or
%                   4. If steps = N, then length(t) == length(w) >= N.
%                   If more values are given in t and w than necessary, the
%                   first N will be used.
%               start - Value of t and w at which to start.
%
%       Outputs:
%               yt - Row vector of symbolic functions in terms of wi+1 to
%                   be solved later.
%                   i.e. wi+1 = f(wi+1)
%

assert(nargin > 5, 'Not enough parameters.');
assert(length(t) == length(w), 'Invalid parameters t and w: Must be equal lengths.');
assert(steps > 1 && steps < 5, 'Invalid step size.');
assert(length(t) >= steps, 'Not enough values given for step size.');

syms w1;

step = ( upperLimit - lowerLimit ) / ( length(t) - 1 );
yt = zeros(1, length(w));

j = steps;
yt(1 : j) = w(1 : j);

for i = j : length(t) - 1
    switch steps
        case 2
            wj = yt(i) + step / 12 * ( 5 * feval(fty, t(i) + step, w1) + 8 * feval(fty, t(i), yt(i)) - feval(fty, t(i - 1), yt(i - 1)));
        case 3
            wj = yt(i) + step / 24 * ( 9 * feval(fty, t(i) + step, w1) + 19 * feval(fty, t(i), yt(i)) - 5 * feval(fty, t(i - 1), yt(i - 1)) + feval(fty, t(i - 2), yt(i - 2)));
        case 4
            wj = yt(i) + step / 720 * ( 251 * feval(fty, t(i) + step, w1) + 646 * feval(fty, t(i), yt(i)) - 264 * feval(fty, t(i - 1), yt(i - 1)) + 106 * feval(fty, t(i - 2), yt(i - 2)) - 19 * feval(fty, t(i - 3), yt(i - 3)));
    end
    yt(i + 1) = vpa(solve(strcat('w1 = ', char(wj))), 15);
    fprintf('y(%5.3f) = %s\n', t(i + 1), yt(i + 1));
end

end
