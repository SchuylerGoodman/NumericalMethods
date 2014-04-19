function [F all] = dividedDiff( xn, yn, type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic

if  strcmpi(type, 'forward')
    forward = 1;
elseif strcmpi(type, 'backward')
    forward = 0;
else
    error('dividedDiff:InvalidInput', 'Input type must be either ''forward'' or ''backward''.  %s was input.', type);
end

len = length(xn);
all = zeros(len, len);
syms x;
F = 0;

for i = 1 : len
    for j = 1 : i
        if j == 1
            fij = yn(i);
        else
            Fij1 = all(i, j - 1);
            Fi1j1 = all(i - 1, j - 1);
            xi = xn(i);
            xij = xn(i - j + 1);
            fij = (Fij1 - Fi1j1) / (xi - xij);
        end
        all(i, j) = fij;
        Fj = fij;
        if forward
            if i == j
                for k = 1 : i - 1
                    Fj = Fj * (x - xn(k));
                end
                F = F + Fj;
            end
        else
            if i == len
                for k = len : -1 : (len - j + 2)
                    Fj = Fj * (x - xn(k));
                end
                F = F + Fj;
            end
        end
    end
end
F = expand(F);
F = matlabFunction(F);
toc
end

