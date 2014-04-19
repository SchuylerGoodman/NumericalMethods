function [res] = findRoot(f, a, b, tol, N)
tic
res = NaN;
if a == b
    disp('Bounds of range cannot be equal.');
    return;
end;

df_a = (f(a + tol) - f(a)) / tol;
df_b = (f(b + tol) - f(b)) / tol;
arr = [a, b];
a = min(arr);
b = max(arr);

if isnan(res) && abs(df_a) >= tol
    res = Newton1(f, a, tol, N);
    disp('Using Newton''s Method');
end
if isnan(res) && abs(df_b) >= tol
    res = Newton1(f, b, tol, N);
    disp('Using Newton''s Method');
end
if isnan(res) && (f(a) <= a && f(b) >= b) || (f(a) >= b && f(b) <= a)
    res = fixedpoint(f, (a + b) / 2, tol, N);
    disp('Using Fixed Point Method');
end
if isnan(res)
    res = bisect_2(f, a, b, tol, N);
    disp('Using Bisection');
end
toc

function [p]=Newton1(f, p0, tol, N)
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

function p=fixedpoint(f, p0, tol, N)
iterations = 0;

while 1
  iterations = iterations + 1;
  if iterations > N
        disp('Cannot find root, may not exist.');
        p = NaN;
        break;
  end
  df = (f(p0 + tol) - f(p0)) / tol;

  if abs(df) < tol
      df = df - tol;
  end
  p = p0 - f(p0) / df;

  if abs(p - p0) < tol
      break; 
  elseif abs(p - p0) > 100 / tol,  
      fprintf('diverges.\n');
      break;
  end; 
  p0 = p;
end
iterations

function p=bisect_2(f,a,b,tol,N)
iterations = 0;
f_a = f(a);
f_b = f(b);
if f_a * f_b >= 0
    disp('Error in bisect_2: f(a) and f(b) have the same sign.');
    p = NaN;
    return;
end
while 1
    iterations = iterations + 1;
    if iterations > N
        disp('Cannot find root, may not exist.');
        p = NaN;
        break;
    end
    p=(a + b) / 2;
    if p - a < tol
        break;
    end
    f_a = f(a);
    f_b = f(b);
    f_p = f(p);
    if f_a * f_p > 0
        a = p; 
    else
        b = p;
    end
end
iterations