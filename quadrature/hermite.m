function [F all] = hermite(a, f, fp)

len = length(a);
all = zeros(2 * len, 2 * len);

syms x;
F = 0;

for i = 1 : 2 * len
    for j = 1 : i
        if j == 1 % If in first column, fill with values from a
            fij = f(round(i / 2));
        elseif j == 2 && ( mod(i, 2) == 0 )
            fij = fp(round(i / 2)); % If in second column, fill with values from f, but only on even i's
        else % Otherwise calculate from previous values
            Fij1 = all(i, j - 1);
            Fi1j1 = all(i - 1, j - 1);
            xi = a(round(i / 2));
            xij = a(round((i - j + 1) / 2));
            fij = (Fij1 - Fi1j1) / (xi - xij);
        end
        all(i, j) = fij; % Fill matrix
        Fj = fij; % Initialize term for this part of the function
        if i == j % Build output function
            for k = 1 : i - 1
                Fj = Fj * (x - a(round(k / 2)));
            end
            F = F + Fj;
        end
    end
end
F = expand(F); % Expand polynomial for easier readability
F = matlabFunction(F); % Convert to function handle for solving

end