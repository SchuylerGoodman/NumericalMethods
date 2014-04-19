function [fdivided, rem] = syndiv(a, root)

% Performs synthetic division on a polynomial.
%
%   input a
%       row vector containing the parent polynomial
%   input root
%       root of a by which syndiv is dividing
%
%   output fdivided
%       row vector containing divided polynomial
%   output rem
%       remainder of the division
%           used for things like Horner's method(xk+1 = xk - P(xk) / Q(xk))

len = length(a);
fdivided = zeros(1, len - 1);
next = 0;

try
    polyval(a, 0);
catch err
    error('MATLAB:syndiv:incorrectInput', ...
        'Input must be a polynomial (row vector).');
end

if len > 1
    for i = 1 : len - 1
        next = a(i) + ( root * next );
        fdivided(i) = next;
    end
end
rem = a(len) + ( root * next );

end

