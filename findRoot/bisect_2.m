
function p=bisect_2(f,a,b,tol,N)
tic
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
toc