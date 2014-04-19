function pxk = lagrange(a, f, xk)
tic
%Uses the Lagrange interpolating polynomial to estimate a value at a given
%x
%   input a
%       row vector containing known values of the function
%   input f
%       either:
%           a function handle representing the function or,
%           a row vector containing the values of f for every x in a
%               length(a) MUST EQUAL length(f)
%   input xk
%       value of f to estimate
%   output pxk
%       estimated value of f(xk)
len = length(a);
pxk = 0;
row = 0;

for i = 1 : len
    li = 1;
    for j = 1 : len
        if i ~= j
            num = xk - a(j);
            den = a(i) - a(j);
            li = li * (num / den);
        end
    end
    try
        fi = feval(f, a(i)); % f must be either a function handle or ...
    catch err
        if (strcmp(err.identifier, 'MATLAB:feval:argMustBeStringOrHandle'))
            try 
                if (row == 0)
                    disp('Input f not a function handle, interpreting as row vector...');
                end
                row = 1;
                fi = f(i); % a row vector of values corresponding to f(a(i)) where i is from 1 to length(a)
            catch err2
                if (strcmp(err2.identifier, 'MATLAB:badsubscript'))
                    msg = sprintf('Attempt to access f(%d) failed; index out of bounds because length(a) ~= length(f).\nLengths must match.\n', i);
                    error('MATLAB:lagrange:mismatch', msg);
                else
                    disp('Found length mismatch.');
                    rethrow(err2);
                end
            end
        else
            rethrow(err);
        end
    end
    
    li = li * fi;
    fprintf('f(x%d) * L%d(%d):    %6.5f\n', i - 1, i - 1, xk, li);
    if i ~= len
        fprintf('                +\n');
    end
    pxk = pxk + li;
end
toc
end

