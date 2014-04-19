function f2x = three_midpoint(x, fx, x0)

%THREE_MIDPOINT  Calculates f'(x0) from given values of x and f(x) using
%                the Three-Point Midpoint Formula for derivatives.
%                1/2h(f(x0+h)-f(x0-h))
%                Needs at least 1 value on each side of x0,
%
%       inputs:
%               x      A vector of EVENLY-SPACED values of x at which f(x)
%                      have been calculated.
%               fx     A vector of values of f(x) corresponding to values
%                      x from the x vector.
%               x0     Value of x at which to calculate f'(x). Must be in x
%                      vector input.
%
%       outputs:
%               f2x    Calculated f'(x0)
%

index = find(x==x0);

if isempty(index)
    error('THREE_MIDPOINT:InvalidInput', ...
        'Value not found in known values.');
end

len = length(x);
lenf = length(fx);

if len < 3 || lenf < 3
    error('THREE_MIDPOINT:InvalidInput', ...
        'It''s called THREE_MIDPOINT because you need at least three values...');
end

if len ~= lenf
    error('THREE_MIDPOINT:InvalidInput', ...
        'First two parameters must be row vectors of equal length.');
end

if abs( len - index ) < 1 || ( index ) <= 1
    error('THREE_MIDPOINT:WrongFormula', ...
        'Not enough values to use midpoint, try THREE_ENDPOINT.');
end

h = x(2) - x(1);
factor = (1 / ( 2 * h ) );

f2x = fx(index + 1) - fx(index - 1);
f2x = factor * f2x;

end
