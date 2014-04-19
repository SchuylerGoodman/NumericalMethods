function [x] = horner(P, x0, tol, iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic
disp('i       x               err');
len = length(P);
i = 0;
while i < iter
    y = P(1);
    z = y;
    x = 0;
    if len > 2
        for nm = 2 : len - 1
            y = x0 * y + P(nm);
            z = x0 * z + y;
        end
    end
    y = x0 * y + P(len);
    if z == 0
        fprintf('Error: division by 0 ==> %.5f - ( %.5f / %.5f )\n', x0, y, z);
        break;
    end
    x = x0 - y / z;
    fprintf('%1d       %6.5f         %6.5f\n', i, x, abs(x - x0));
    if abs(x - x0) < tol
        break;
    end
    i = i + 1;
    x0 = x;
end
toc
end

