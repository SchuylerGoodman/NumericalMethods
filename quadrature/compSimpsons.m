function [app] = compSimpsons(f, left, right, n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if mod(n, 2) ~= 0
    error('Input n must be an even number: %d given', n);
end

step = ( right - left ) / n;

fA = f(left);
fB = f(right);

f1 = 0;
f2 = 0;

for i = 1 : n - 1
    x = left + i * step;
    if mod(i, 2) == 1
        f2 = f2 + feval(f, x);
    else
        f1 = f1 + feval(f, x);
    end
end

app = ( step / 3 ) * ( fA + 2 * f1 + 4 * f2 + fB);

end

