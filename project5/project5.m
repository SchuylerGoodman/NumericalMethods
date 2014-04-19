function [] = project5()
%PROJECT5 Code for Project 5 from Math 410 - Intro to Numerical Methods.
%   Assignment prompt included in the folder.

U1axp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U1bx;
U1bxp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U1ax-U2ax)^2/(sqrt((U1ax-U2ax)^2+(U1ay-U2ay)^2))^3 + (U1ax-U3ax)^2/(sqrt((U1ax-U3ax)^2+(U1ay-U3ay)^2))^3 + (U1ax-U4ax)^2/(sqrt((U1ax-U4ax)^2+(U1ay-U4ay)^2))^3;
U1ayp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U1by;
U1byp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U1ay-U2ay)^2/(sqrt((U1ax-U2ax)^2+(U1ay-U2ay)^2))^3 + (U1ay-U3ay)^2/(sqrt((U1ax-U3ax)^2+(U1ay-U3ay)^2))^3 + (U1ay-U4ay)^2/(sqrt((U1ax-U4ax)^2+(U1ay-U4ay)^2))^3;

U2axp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U2bx;
U2bxp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U2ax-U1ax)^2/(sqrt((U2ax-U1ax)^2+(U2ay-U1ay)^2))^3 + (U2ax-U3ax)^2/(sqrt((U2ax-U3ax)^2+(U2ay-U3ay)^2))^3 + (U2ax-U4ax)^2/(sqrt((U2ax-U4ax)^2+(U2ay-U4ay)^2))^3;
U2ayp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U2by;
U2byp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U2ay-U1ay)^2/(sqrt((U2ax-U1ax)^2+(U2ay-U1ay)^2))^3 + (U2ay-U3ay)^2/(sqrt((U2ax-U3ax)^2+(U2ay-U3ay)^2))^3 + (U2ay-U4ay)^2/(sqrt((U2ax-U4ax)^2+(U2ay-U4ay)^2))^3;

U3axp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U3bx;
U3bxp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U3ax-U1ax)^2/(sqrt((U3ax-U1ax)^2+(U3ay-U1ay)^2))^3 + (U3ax-U2ax)^2/(sqrt((U3ax-U2ax)^2+(U3ay-U2ay)^2))^3 + (U3ax-U4ax)^2/(sqrt((U3ax-U4ax)^2+(U3ay-U4ay)^2))^3;
U3ayp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U3by;
U3byp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U3ay-U1ay)^2/(sqrt((U3ax-U1ax)^2+(U3ay-U1ay)^2))^3 + (U3ay-U2ay)^2/(sqrt((U3ax-U2ax)^2+(U3ay-U2ay)^2))^3 + (U3ay-U4ay)^2/(sqrt((U3ax-U4ax)^2+(U3ay-U4ay)^2))^3;

U4axp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U4bx;
U4bxp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U4ax-U1ax)^2/(sqrt((U4ax-U1ax)^2+(U4ay-U1ay)^2))^3 + (U4ax-U2ax)^2/(sqrt((U4ax-U2ax)^2+(U4ay-U2ay)^2))^3 + (U4ax-U3ax)^2/(sqrt((U4ax-U3ax)^2+(U4ay-U3ay)^2))^3;
U4ayp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) U4by;
U4byp = @(t, U1ax, U1ay, U2ax, U2ay, U3ax, U3ay, U4ax, U4ay, U1bx, U1by, U2bx, U2by, U3bx, U3by, U4bx, U4by) (U4ay-U1ay)^2/(sqrt((U4ax-U1ax)^2+(U4ay-U1ay)^2))^3 + (U4ay-U2ay)^2/(sqrt((U4ax-U2ax)^2+(U4ay-U2ay)^2))^3 + (U4ay-U3ay)^2/(sqrt((U4ax-U3ax)^2+(U4ay-U3ay)^2))^3;


fty = {U1axp, U1bxp, U1ayp, U1byp, U2axp, U2bxp, U2ayp, U2byp, U3axp, U3bxp, U3ayp, U3byp, U4axp, U4bxp, U4ayp, U4byp};
lower = 0;
upper = 3;
initials = [1.0597 1.7696 0 -.8094 -1.0597 1.7696 0 -2.7299 -.5537 -.3988 1.0934 0.0000 -.5539 .3989 .0142 0.000];
intervals = 100;

[t, wt] = rk4sys(fty, lower, upper, initials, intervals);

figure
plot(wt(1, :), wt(2, :), '-r', wt(3, :), wt(4, :), '-b', wt(5, :), wt(6, :), '-g', wt(7, :), wt(8, :), '-m');
legend('q1(t)', 'q2(t)', 'q3(t)', 'q4(t)');

end

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
%fprintf('tj        ');
% for k = 1 : numEqs
%     fprintf('w%dj            ', k);
% end
% fprintf('\n');

%  Initialize return values
t = zeros(1, numPoints);
t(1) = lowerLimit;
wt = zeros(numEqs, numPoints);
% fprintf('%3.2f      ', t(1));  %  Display first row
for init = 1 : numEqs
    wt(init, 1) = initials(init);
%     fprintf('%6.7f      ', wt(init, 1));
end
% fprintf('\n');

%  Initialize intermediate values
K1 = zeros(1, numEqs);
K2 = zeros(1, numEqs);
K3 = zeros(1, numEqs);
K4 = zeros(1, numEqs);

%  Iterate through mesh points and calculate estimates with rk4
for i = 2 : numPoints
    t(i) = t(i - 1) + step;  %  Increment t
%     fprintf('%3.2f      ', t(i));  %  Print t
    
    %  Confusing notation for calculations because K1, K2, K3, and K4 are
    %  row vectors, but the required values in wt are in a column.
    K1wt = wt(:, i - 1)';
    for j = 1 : numEqs
        K1(j) = step * f(fty, j, t(i - 1), K1wt);
    end
    
    K2wt = wt(:, i - 1)' + 1 / 2 * K1;
    for j = 1 : numEqs
        K2(j) = step * f(fty, j, t(i - 1) + step / 2, K2wt);
    end
    
    K3wt = wt(:, i - 1)' + 1 / 2 * K2;
    for j = 1 : numEqs
        K3(j) = step * f(fty, j, t(i - 1) + step / 2, K3wt);
    end
    
    K4wt = K1wt + K3;
    for j = 1 : numEqs
        K4(j) = step * f(fty, j, t(i - 1) + step, K4wt);
    end
    
    for j = 1 : numEqs
        wt(j, i) = wt(j, i - 1) + ( K1(j) + 2 * K2(j) + 2 * K3(j) + K4(j) ) / 6;
%         fprintf('%6.7f      ', wt(j, i));  %  Print estimates
    end
%     fprintf('\n');
    
end

end

function [Kj] = f(fty, j, t, w)
%F          Helper function for rk4sys.  Calculates value of the function
%               at j with inputs t and all values in w.  Assumes that given
%               function at fty{j} needs length(w) + 1 inputs.

    wCell = num2cell(w);
    Kj = feval(fty{j}, t, wCell{:});

end