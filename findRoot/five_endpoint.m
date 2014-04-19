function f2x = five_endpoint(x, fx, x0)

%FIVE_ENDPOINT   Calculates f'(x0) from given values of x and f(x) using
%                the Five-Point Endpoint Formula for derivatives.
%                1/12h(-25f(x0)+48f(x0+h)-36f(x0+2h)+16f(x0+3h)-3f(x0+4h))
%                Automatically determines whether to use Forward-difference
%                or Backward-difference.  Needs at least five values before
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
    error('FIVE_ENDPOINT:InvalidInput', ...
        'Value not found in known values.');
end

len = length(x);
lenf = length(fx);

if len < 5 || lenf < 5
    error('FIVE_ENDPOINT:InvalidInput', ...
        'It''s called FIVE_ENDPOINT because you need at least five values...');
end

if len ~= lenf
    error('FIVE_ENDPOINT:InvalidInput', ...
        'First two parameters must be row vectors of equal length.');
end

if ( len - index ) < 4 && ( index ) <= 4
    error('FIVE_ENDPOINT:WrongFormula', ...
        'Not enough values to use endpoint, try FIVE_MIDPOINT.');
end

h = x(2) - x(1);
factor = (1 / ( 12 * h ) );

if ( len - index ) >= 4
    f2x = -25 * fx(index) + 48 * fx(index + 1) - 36 * fx(index + 2) ...
        + 16 * fx(index + 3) - 3 * fx(index + 4);
elseif index > 4
    f2x = -25 * fx(index) + 48 * fx(index - 1) - 36 * fx(index - 2) ...
        + 16 * fx(index - 3) - 3 * fx(index - 4);
else
    error('FIVE_ENDPOINT:UselessError', ...
        'I have no idea what went wrong, because I''m too tired to think about it.');
end
f2x = factor * f2x;

end

