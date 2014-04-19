function p=fixedpoint(f, p0, tol, N)
tic
iterations = 0;

%syms x;
%g = sym(f);
%df = diff(g, x, 1);
%f = x - g / df;

while 1
  iterations = iterations + 1;
  if iterations > N
        disp('Cannot find root, may not exist.');
        p = NaN;
        break;
  end
  %dfp = subs(df, x, p0);
  df = (f(p0 + tol) - f(p0)) / tol;

  %if abs(dfp) < tol
  if abs(df) < tol
      %dfp = dfp - tol;
      df = df - tol;
  end
  %p = p0 - subs(f, x, p0) / dfp;
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
toc