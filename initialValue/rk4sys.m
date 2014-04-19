function [t, wt] = rk4sys(fty, lowerLimit, upperLimit, initials, intervals)
%RK4SYS Estimate value of a system of initial value problems at the given
%               upper limit using Runge-Kutta Method of order 4.
%
%       Inputs:
%               fty - f(t, y), or y'.  The cell array containing the system
%                   of differential equations whose solutions are being
%                   estimated.
%               lowerLimit - Lower bound for t
%               upperLimit - Upper bound for t, value at which to estimate
%                   y(t)
%               initials - Row vector of initial y(lowerLimit) = initials
%                   for all equations in fty.
%               intervals - Number of desired mesh points - 1
%
%       Outputs:
%               t - Matrix containing row vectors of mesh points used for
%                   each equation in fty.
%               wt - Matrix containing row vectors of estimated values of
%                   y(t) at mesh points for each equation in fty.
%               wdt - Matrix containing row vectors of values of y'(t)
%                   using estimated values of y(t) stored in yt for each
%                   equation in fty.
%               NOTE:
%                   For all outputs, each row corresponds to an equation
%                   from the system in fty and each column corresponds to a
%                   mesh point.
%
%       NOTE:
%               t, yt, and ydt can be easily used with HERMITE, CUBIC_NAT, and
%                   CUBIC_CLAMP to interpolate other values of y(t).
%

assert(isa(fty, 'cell'), 'Given system not a cell array.');
assert(isa(lowerLimit, 'double'), 'Given lower limit not a number.');
assert(length(lowerLimit) == 1, 'Given lower limit is an array.');
assert(isa(upperLimit, 'double'), 'Given upper limit not a number.');
assert(length(upperLimit) == 1, 'Given upper limit is an array.');
assert(lowerLimit < upperLimit, 'Given limits out of order.');
assert(isa(initials, 'double'), 'Given initial values array contains non-numbers.');
assert(length(fty) <= length(initials), 'Not enough initial values for given system.');
assert(isa(intervals, 'double'), 'Given number of intervals contains non-number.');
assert(length(intervals) == 1, 'Given number of intervals is an array.');

step = ( upperLimit - lowerLimit ) / intervals;
numPoints = intervals + 1;
numEqs = length(fty);

%  Display table header
fprintf('tj        ');
for k = 1 : numEqs
    fprintf('w%dj            ', k);
end
fprintf('\n');

%  Initialize return values
t = zeros(1, numPoints);
t(1) = lowerLimit;
wt = zeros(numEqs, numPoints);
fprintf('%3.2f      ', t(1));  %  Display first row
for init = 1 : numEqs
    wt(init, 1) = initials(init);
    fprintf('%6.7f      ', wt(init, 1));
end
fprintf('\n');

%  Initialize intermediate values
K1 = zeros(1, numEqs);
K2 = zeros(1, numEqs);
K3 = zeros(1, numEqs);
K4 = zeros(1, numEqs);

%  Iterate through mesh points and calculate estimates with rk4
for i = 2 : numPoints
    t(i) = t(i - 1) + step;  %  Increment t
    fprintf('%3.2f      ', t(i));  %  Print t
    
    %  Confusing notation for calculations because K1, K2, K3, and K4 are
    %  row vectors, but the required values in wt are in a column.
    for j = 1 : numEqs
        K1(j) = step * f(fty, j, t(i - 1), wt(:, i - 1)');
    end
    
    for j = 1 : numEqs
        K2(j) = step * f(fty, j, t(i - 1) + step / 2, wt(:, i - 1)' + 1 / 2 * K1);
    end
    
    for j = 1 : numEqs
        K3(j) = step * f(fty, j, t(i - 1) + step / 2, wt(:, i - 1)' + 1 / 2 * K2);
    end
    
    for j = 1 : numEqs
        K4(j) = step * f(fty, j, t(i - 1) + step, wt(:, i - 1)' + K3);
    end
    
    for j = 1 : numEqs
        wt(j, i) = wt(j, i - 1) + ( K1(j) + 2 * K2(j) + 2 * K3(j) + K4(j) ) / 6;
        fprintf('%6.7f      ', wt(j, i));  %  Print estimates
    end
    fprintf('\n');
    
end

end

function [Kj] = f(fty, j, t, w)
%F          Helper function for rk4sys.  Calculates value of the function
%               at j with inputs t and w{:}.

    wCell = num2cell(w);
    Kj = feval(fty{j}, t, wCell{:});

end