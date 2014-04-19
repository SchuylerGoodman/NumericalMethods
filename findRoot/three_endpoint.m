function f2x = three_endpoint(x, fx, x0)

%THREE_ENDPOINT  Calculates f'(x0) from given values of x and f(x) using
%                the Three-Point Endpoint Formula for derivatives.
%                1/2h(-3f(x0)+4f(x0+h)-f(x0+2h))
%                Automatically determines whether to use Forward-difference
%                or Backward-difference.  Needs at least three values before
%                or after the given x0.
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
    error('THREE_ENDPOINT:InvalidInput', ...
        'Value not found in known values.');
end

len = length(x);
lenf = length(fx);

if len < 3 || lenf < 3
    error('THREE_ENDPOINT:InvalidInput', ...
        'It''s called THREE_ENDPOINT because you need at least three values...');
end

if len ~= lenf
    error('THREE_ENDPOINT:InvalidInput', ...
        'First two parameters must be row vectors of equal length.');
end

if ( len - index ) < 2 && ( index ) <= 2
    error('THREE_ENDPOINT:WrongFormula', ...
        'Not enough values to use endpoint, try THREE_MIDPOINT.');
end

h = x(2) - x(1);
factor = (1 / ( 2 * h ) );

if ( len - index ) >= 2
    f2x = -3 * fx(index) + 4 * fx(index + 1) - fx(index + 2);
elseif index > 2
    f2x = -3 * fx(index) + 4 * fx(index - 1) - fx(index - 2);
else
    error('THREE_ENDPOINT:UselessError', ...
        'I have no idea what went wrong, because I''m too tired to think about it.');
end
f2x = factor * f2x;

end

