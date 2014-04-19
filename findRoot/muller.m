function p = muller(f, p0, p1, p2, tol, iter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tic
i = 3;
E = 0;
disp('i       x               err');
while i < iter
    h1 = p1 - p0;
    h2 = p2 - p1;
    d1 = (polyval(f, p1) - polyval(f, p0)) / h1;
    d2 = (polyval(f, p2) - polyval(f, p1)) / h2;
    d = (d2 - d1) / (h2 + h1);
    b = d2 + h2 * d;
    D = sqrt(b^2 - 4 * polyval(f, p2) * d);
    
    if abs(b - D) < abs(b + D)
        E = b + D;
    else
        E = b - D;
    end
    h = -2 * polyval(f, p2) / E;
    p = p2 + h;
    fprintf('%1d       %6.5f         %6.5f\n', i, p, abs(h));
    if abs(h) < tol
        break;
    end
    p0 = p1;
    p1 = p2;
    p2 = p;
    i = i + 1;
end
toc
end

