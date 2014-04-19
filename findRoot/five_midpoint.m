function f2x = five_midpoint(x, fx, x0)

%FIVE_MIDPOINT   Calculates f'(x0) from given values of x and f(x) using
%                the Five-Point Midpoint Formula for derivatives.
%                1/12h(-25f(x0)+48f(x0+h)-36f(x0+2h)+16f(x0+3h)-3f(x0+4h))
%                Needs at least 2 values on each side of x0,
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
    error('FIVE_MIDPOINT:InvalidInput', ...
        'Value not found in known values.');
end

len = length(x);
lenf = length(fx);

if len < 5 || lenf < 5
    error('FIVE_MIDPOINT:InvalidInput', ...
        'It''s called FIVE_MIDPOINT because you need at least five values...');
end

if len ~= lenf
    error('FIVE_MIDPOINT:InvalidInput', ...
        'First two parameters must be row vectors of equal length.');
end

if abs( len - index ) < 2 || ( index ) <= 2
    error('FIVE_MIDPOINT:WrongFormula', ...
        'Not enough values to use midpoint, try FIVE_ENDPOINT.');
end

h = x(2) - x(1);
factor = (1 / ( 12 * h ) );

f2x = fx(index - 2) - 8 * fx(index - 1) + 8 * fx(index + 1) - fx(index + 2);
f2x = factor * f2x;

end
