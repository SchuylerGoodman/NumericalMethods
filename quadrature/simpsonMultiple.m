function value = simpsonMultiple(f, lower, upper, intervals, multiplicity)
%SIMPSONMULTIPLE Only works for up to triple integrals
tic
level = 0;
value = simpsonMultipleRecursive(f, lower, upper, intervals, multiplicity, level + 1);
toc
end

function value = simpsonMultipleRecursive(f, lower, upper, intervals, multiplicity, level)

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
if (level > multiplicity)
    throw('SIMPSONRECURSIVE', 'On level %d when multiplicity given is %d.', level, multiplicity);
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
    step = ( feval(upper{level}, 0) - feval(lower{level}, 0) ) / intervals(level);
    nextLow = feval(lower{level}, 0);
    nextUp = feval(upper{level}, 0);
elseif level == 2
    step = ( feval(upper{level}, lower{1}) - feval(lower{level}, lower{1}) ) / intervals(level);
    nextLow = feval(lower{level}, lower{1});
    nextUp = feval(upper{level}, upper{1});
else
    step = ( feval(upper{level}, lower{1}, lower{2}) - feval(lower{level}, lower{1}, lower{2}) ) / intervals(level);
    nextLow = feval(lower{level}, lower{1}, lower{2});
    nextUp = feval(upper{level}, upper{1}, upper{2});
end

ends = 0;
evens = 0;
odds = 0;

for count = 1 : intervals(level) + 1
    tempLower{level} = nextLow + ( count - 1 ) * step;
    tempUpper{level} = nextUp + ( count - 1 ) * step;
    if level < multiplicity
        recValue = simpsonMultipleRecursive(f, tempLower, tempUpper, intervals, multiplicity, level + 1);
    else
        if level == 1
            recValue = feval(f, tempLower{level});
        elseif level == 2
            recValue = feval(f, tempLower{level - 1}, tempLower{level});
        else
            recValue = feval(f, tempLower{level - 2}, tempLower{level - 1}, tempLower{level});
        end
    end
        
    if count - 1 == 0 || count - 1 == intervals(level)
        ends = ends + recValue;
    elseif mod(count - 1, 2) % if count - 1 is odd
        odds = odds + recValue;
    else % if count - 1 is even
        evens = evens + recValue;
    end
end

value = ( step / 3 ) * ( ends + 2 * evens + 4 * odds );

end