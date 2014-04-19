function y = solve_spline( cspl, x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic
[n m] = size(cspl);
npts = length(x);
y = zeros(1, npts);

for i = 1 : npts
    row = max( find( cspl(1 : n - 1) < x(i) ) );
    if isempty(row)
        row = 1;
    end
    
    xdiff = x(i) - cspl(row);
    temp = cspl(row, 2);
    temp = temp + cspl(row, 3) * xdiff;
    if m == 5
        temp = temp + ( cspl(row, 4) * xdiff ^ 2 );
        temp = temp + ( cspl(row, 5) * xdiff ^ 3 );
    end
    y(1, i) = temp;
end
toc
end

