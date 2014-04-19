function [p]=Newton1(f, p0, tol, N)
tic
iterations = 0;
p = NaN;
while 1
    iterations = iterations + 1;
    if iterations > N
        disp('Cannot find root, may not exist.');
        p = NaN;
        break;
    end
    df = (f(p0 + tol) - f(p0)) / tol;
    count = 0;
    while abs(df) == 0 % if df is 0
        if count > N
	        disp('Newton''s Method failed: df is 0');
            toc
            return;
        end
        if abs(f(p0)) >= tol % and p0 is not the root
            p0 = p0 - tol * 10^count; % move p0 tol to the left
        else
            break;
        end
        df = (f(p0 + tol) - f(p0)) / tol; % reevaluate df and try again
        count = count + 1;
    end
p = p0 - f(p0) / df;
if abs(p - p0) < tol
    break;
end
p0=p;
end
iterations
toc