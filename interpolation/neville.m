function [all, x] = neville(a, f, xk)
%Uses Neville's method to estimate a value at a given x
%   input a
%       row vector containing known values of the function
%   input f
%       either:
%           a function handle representing the function or,
%           a row vector containing the values of f for every x in a
%               length(a) MUST EQUAL length(f)
%   input xk
%       value of f to estimate
%   output all
%       matrix with all recursive results
%   output x
%       estimated value of f(xk)
len = length(a);
all = zeros(len, len);
row_vec = 0;

for i = 1 : len
    for j = 1 : i
        fij = 0;
        if (j == 1)
            try
                fij = feval(f, a(i));
            catch err
                if (strcmp(err.identifier, 'MATLAB:feval:argMustBeStringOrHandle'))
                    try 
                        if (row_vec == 0)
                            disp('Input f not a function handle, interpreting as row vector...');
                        end
                        row_vec = 1;
                        fij = f(i); % a row vector of values corresponding to f(a(i)) where i is from 1 to length(a)
                    catch err2
                        if (strcmp(err2.identifier, 'MATLAB:badsubscript'))
                            msg = sprintf('Attempt to access f(%d) failed; index out of bounds because length(a) ~= length(f).\nLengths must match.\n', i);
                            error('MATLAB:lagrange:mismatch', ...
                                'Attempt to access f(%d) failed; index out of bounds because length(a) ~= length(f).\nLengths must match.\n', i ...
                            );
                        else
                            disp('Found length mismatch.');
                            rethrow(err2);
                        end
                    end
                else
                    rethrow(err);
                end
            end
        else
            xd0 = xk - a(i - j + 1);
            Q1 = all(i, j - 1);
            xd1 = xk - a(i);
            Q0 = all(i - 1, j - 1);
            xdb = a(i) - a(i - j + 1);
            fij = ( ( xd0 * Q1 ) - ( xd1 * Q0 ) ) / ( xdb );
        end
        all(i, j) = fij;
    end
    x = all(i, i);
end
        

end

