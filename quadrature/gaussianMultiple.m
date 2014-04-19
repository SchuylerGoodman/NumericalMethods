function value = gaussianMultiple(f, lower, upper, intervals, multiplicity)
%GAUSSIANMULTIPLE   Only works for up to triple integrals
%                   The input parameters for this function are more than
%                   likely needlessly complex.
%   

r = zeros(5, 5);
r(1, :) = 0;
r(2, :) = [1/sqrt(3) -1/sqrt(3) 0 0 0];
r(3, :) = [sqrt(3/5) 0 -sqrt(3/5) 0 0];
r(4, :) = [0.861136311594053 0.339981043584856 -0.339981043584857 -0.861136311594054 0];
r(5, :) = [.9061798459 .5384693101 0 -.5384693101 -.9061798459];

c = zeros(5, 5);
c(1, :) = 0;
c(2, :) = [1 1 0 0 0];
c(3, :) = [5/9 8/9 5/9 0 0];
c(4, :) = [.3478548451 .6521451549 .6521451549 .3478548451 0];
c(5, :) = [.2369268850 .4786286705 .568888888888888 .4786286705 .2369268850];

level = 0;
value = gaussianMultipleRecursive(f, r, c, lower, upper, intervals, multiplicity, level + 1);

end

function value = gaussianMultipleRecursive(f, r, c, lower, upper, intervals, multiplicity, level)

assert(isa(f, 'function_handle'), 'Input must be a of type function handle.');
assert(isa(lower, 'cell'), 'Lower bounds must be in a cell array.');
assert(isa(upper, 'cell'), 'Upper bounds must be in a cell array.');
assert(isa(intervals, 'double'), 'Interval counts must be in a double array.');
assert(isa(multiplicity, 'double') && multiplicity > 0 && multiplicity < 4, ...
       'Multiplicity must be a positive integer less than 4.');
assert(length(lower) == length(upper) && ...
       length(upper) == length(intervals) && ...
       length(intervals) == multiplicity, ...
       'All arrays must have the same length as the multiplicity.');
assert(max(intervals) <= 5, 'Interval cannot be greater than 5, sorry.');
if (level > multiplicity)
    throw('GAUSSIANRECURSIVE', 'On level %d when multiplicity given is %d.', level, multiplicity);
end

assert(isa(lower{level}, 'function_handle') || ...
       isa(lower{level}, 'double'), ...
       'Lower bound must be of type function_handle or double.');
assert(isa(upper{level}, 'function_handle') || ...
       isa(upper{level}, 'double'), ...
       'Upper bound must be of type function_handle or double.');
   
if isa(lower{level}, 'double')
    if level < 3
        lower{level} = @(x) lower{level};
    else
        lower{level} = @(x, y) lower{level};
    end
end
if isa(upper{level}, 'double')
    if level < 3
        upper{level} = @(x) upper{level};
    else
        upper{level} = @(x, y) upper{level};
    end
end

tempLower = lower;
tempUpper = upper;
if level == 1
    nextLow = feval(lower{level}, 0);
    nextUp = feval(upper{level}, 0);
elseif level == 2
    x = lower{level - 1};
    nextLow = feval(lower{level}, x);
    nextUp = feval(upper{level}, x);
else
    x = lower{level - 2};
    y = lower{level - 1};
    nextLow = feval(lower{level}, x, y);
    nextUp = feval(upper{level}, x, y);
end

h1 = ( nextUp - nextLow ) / 2;
h2 = ( nextUp + nextLow ) / 2;
recValue = 0;

for count = 1 : intervals(level)
    tempLower{level} = h1 * r(intervals(level), count) + h2;
    if level < multiplicity
        recValue = recValue + c(intervals(level), count) * gaussianMultipleRecursive(f, r, c, tempLower, tempUpper, intervals, multiplicity, level + 1);
    else
        if level == 1
            recValue = recValue + c(intervals(level), count) * feval(f, tempLower{level});
        elseif level == 2
            recValue = recValue + c(intervals(level), count) * feval(f, tempLower{level - 1}, tempLower{level});
        else
            recValue = recValue + c(intervals(level), count) * feval(f, tempLower{level - 2}, tempLower{level - 1}, tempLower{level});
        end
    end
end

value = h1 * recValue;

end