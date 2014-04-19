function [] = project4(xi)

%PROJECT4  Compute the estimate of input xi with varying degrees of
%          accuracy using both natural and clamped cubic splines, as well
%          as an calculating the hermite polynomial and an estimate at
%          x=2.8.
%
%     calling sequences:
%             project4(xi);
%
%     inputs:
%             xi      Numeric value between 0 and 3
%
%     outputs:
%             NONE
%
%     NOTES:
%             No output is returned, but a large amount of text is printed,
%             as well as several plots.

if xi < 0 || xi > 3
    error('PROJECT4:InvalidParameter', 'Must input a value between 0 and 3: %6.5f was input.', xi);
end

tic
syms x;
f = @(x) exp(x^2 / 2) * sin(x);

df = matlabFunction(diff(f, x));
trueAns = f(xi);

spline0 = 0 : .5 : 3;
fspline0 = zeros(1, length(spline0));
spline1 = 0 : .3 : 3;
fspline1 = zeros(1, length(spline1));
spline2 = 0 : .25 : 3;
fspline2 = zeros(1, length(spline2));
spline3 = [0 .5 1 1.5:.2:2.3 2.4:.1:3];
fspline3 = zeros(1, length(spline3));

for i = 1 : length(spline0)
    fspline0(1, i) = f(spline0(1, i));
end
for i = 1 : length(spline1)
    fspline1(1, i) = f(spline1(1, i));
end
for i = 1 : length(spline2)
    fspline2(1, i) = f(spline2(1, i));
end
for i = 1 : length(spline3)
    fspline3(1, i) = f(spline3(1, i));
end
fO = df(0);
fN = df(3);

%===========NATURAL SPLINES================================================
% disp('NATURAL SPLINES');
% fprintf('Worst even natural spline: %s\n', num2str(spline0));
% disp('X         f(X)              Error');
sol0 = solve_spline(cubic_nat(spline0, fspline0), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, sol0, abs(sol0 - trueAns));

% fprintf('Second best even natural spline: %s\n', num2str(spline1));
% disp('X         f(X)              Error');
sol1 = solve_spline(cubic_nat(spline1, fspline1), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, sol1, abs(sol1 - trueAns));

% fprintf('Best even natural spline: %s\n', num2str(spline2));
% disp('X         f(X)              Error');
sol2 = solve_spline(cubic_nat(spline2, fspline2), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, sol2, abs(sol2 - trueAns));

% fprintf('Accurate endpoint natural spline: %s\n', num2str(spline3));
% disp('X         f(X)              Error');
sol3 = solve_spline(cubic_nat(spline3, fspline3), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, sol3, abs(sol3 - trueAns));

%===========CLAMPED SPLINES================================================
% disp('CLAMPED SPLINES');
% fprintf('Worst even clamped spline: %s\n', num2str(spline0));
% disp('X         f(X)              Error');
csol0 = solve_spline(cubic_clamp(spline0, fspline0, fO, fN), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, csol0, abs(csol0 - trueAns));

% fprintf('Second best even clamped spline: %s\n', num2str(spline1));
% disp('X         f(X)              Error');
csol1 = solve_spline(cubic_clamp(spline1, fspline1, fO, fN), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, csol1, abs(csol1 - trueAns));

% fprintf('Best even clamped spline: %s\n', num2str(spline2));
% disp('X         f(X)              Error');
csol2 = solve_spline(cubic_clamp(spline2, fspline2, fO, fN), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, csol2, abs(csol2 - trueAns));

% fprintf('Accurate endpoint clamped spline: %s\n', num2str(spline3));
% disp('X         f(X)              Error');
csol3 = solve_spline(cubic_clamp(spline3, fspline3, fO, fN), xi);
% fprintf('%3.2f      %6.5f         %6.5f\n\n', xi, csol3, abs(csol3 - trueAns));

disp('Spline Method  X Values                     X          f(X)         Error');
fprintf('Actual Value   N/A                          %3.2f       %6.5f      %6.5f\n', xi, trueAns, 0);
fprintf('Natural        [0 : .5 : 3]                 %3.2f       %6.5f      %6.5f\n', xi, sol0, abs(sol0 - trueAns));
fprintf('Natural        [0 : .3 : 3]                 %3.2f       %6.5f      %6.5f\n', xi, sol1, abs(sol1 - trueAns));
fprintf('Natural        [0 : .25 : 3]                %3.2f       %6.5f      %6.5f\n', xi, sol2, abs(sol2 - trueAns));
fprintf('Natural        [0 .5 1 1.5:.2:2.3 2.4:.1:3] %3.2f       %6.5f      %6.5f\n', xi, sol3, abs(sol3 - trueAns));
fprintf('Clamped        [0 : .5 : 3]                 %3.2f       %6.5f      %6.5f\n', xi, csol0, abs(csol0 - trueAns));
fprintf('Clamped        [0 : .3 : 3]                 %3.2f       %6.5f      %6.5f\n', xi, csol1, abs(csol1 - trueAns));
fprintf('Clamped        [0 : .25 : 3]                %3.2f       %6.5f      %6.5f\n', xi, csol2, abs(csol2 - trueAns));
fprintf('Clamped        [0 .5 1 1.5:.2:2.3 2.4:.1:3] %3.2f       %6.5f      %6.5f\n\n', xi, csol3, abs(csol3 - trueAns));

xplot = 0 : .001 : 3;
xspline = [0 1 1.6 2.3 2.6 2.8 3];
yspline = exp(xspline.^2 / 2).* sin(xspline);
csplCla = cubic_clamp(xspline, yspline, fO, fN);
yplotReal = exp(xplot.^2 / 2).* sin(xplot);
yplotCla = solve_spline(csplCla, xplot);

disp('Using [0 1 1.6 2.3 2.6 2.8 3] for spline nodes (7 points):');
[m i] = max(abs(yplotReal - yplotCla));
fprintf('Maximum error is %6.5f at x = %4.3f.\n\n', m, i / 1000);
% for i = 1 : length(xplot)
%     yplotNat(1, i) = solve_spline(csplNat, xplot(1, i));
%     yplotCla(1, i) = solve_spline(csplCla, xplot(1, i));
%     yplotReal(1, i) = f(xplot(1, i));
% end

plot(xplot, yplotReal, 'blue');
hold all;
plot(xplot, yplotCla, 'green');
hold off;
% fprintf('Splines plotted...\n\n');

%==========HERMITE=========================================================
dxs = zeros(1, length(spline0));
for i = 1 : length(spline0)
    dxs(1, i) = df(spline0(1, i));
end
[F hermTable] = hermite(spline0, fspline0, dxs);
disp('Hermite Table');
disp(hermTable);
disp('Hermite Polynomial');
disp(F);
disp('Evaluated at x = 2.8: ');
disp(F(2.8));
fprintf('Error is %6.5f\n', abs(F(2.8) - f(2.8)));

toc
end

function cspl = cubic_nat(x, a)

%CUBIC_NAT  compute the cubic spline interpolant, subject to 
%              "not-a-knot" boundary conditions, associated with a 
%              given set of interpolating points and function values
%
%     calling sequences:
%             csn = cubic_nat ( x, a )
%             cubic_nat ( x, a )
%
%     inputs:
%             x       vector containing the interpolating points
%                     (must be sorted in ascending order)
%             a       vector containing function values
%                     the i-th entry in this vector is the function
%                     value associated with the i-th entry in the 'x'
%                     vector
%
%     output:
%             cspl    five column matrix containing the information which
%                     defines the "not-a-knot" cubic spline interpolant
%                     - first column: interpolating points
%                     - second column: function values
%                     - third column: coefficients of linear terms
%                     - fourth column: coefficients of quadratic terms
%                     - fifth column: coefficients of cubic terms
%
%     NOTE:
%             to evaluate the "not-a-knot" cubic spline interpolant apply
%             the routine SOLVE_SPLINE to the output of this routine 
%

n = length(x) - 1;
if length(a) ~= length(x)
    error('cubic_nat:InvalidParameters', 'Length of parameters must match');
end

hdiff = zeros(1, n);
alf = zeros(1, n - 1);
l = zeros(1, n + 1);
m = zeros(1, n);
z = l;
c = z;
b = c;
d = c;

l(1) = 1;
l(n + 1) = 1;
m(1) = 0;
z(1) = 0;
z(n + 1) = 0;
b(n + 1) = NaN;
c(n + 1) = 0;
d(n + 1) = NaN;

for i = 1 : n
    hdiff(i) = x(i + 1) - x(i);
    if i > 1
        alf1 = (3.0 / hdiff(i)) * (a(i + 1) - a(i));
        alf2 = ( 3.0 / hdiff(i - 1)) * ( a(i) - a(i - 1) );
        alf(i - 1) = alf1 - alf2;
        
        l(i) = 2.0 * ( x(i + 1) - x(i - 1)) - ( hdiff(i - 1) * m(i - 1) );
        m(i) = hdiff(i) / l(i);
        z(i) = ( alf(i - 1) - ( hdiff(i - 1) * z(i - 1) ) ) / l(i);
    end
end
for j = n : -1 : 1
    c(j) = z(j) - (m(j) * c(j + 1));
    b(j) = ( a(j + 1) - a(j) ) / hdiff(j) - ( hdiff(j) * ( c(j + 1) + 2.0 * c(j) ) ) / 3.0;
    d(j) = ( c(j + 1) - c(j)) / ( 3.0 * hdiff(j) );
end

if ( nargout == 0 )
    disp( [x' a' b' c' d'] );
else
    cspl = [x' a' b' c' d'];
end

end

function cspl = cubic_clamp( x, fx, fpO, fpN )

%CUBIC_CLAMP  compute the cubic spline interpolant, subject to 
%              "not-a-knot" boundary conditions, associated with a 
%              given set of interpolating points and function values
%
%     calling sequences:
%             csn = cubic_clamp ( x, a )
%             cubic_clamp ( x, a )
%
%     inputs:
%             x       vector containing the interpolating points
%                     (must be sorted in ascending order)
%             a       vector containing function values
%                     the i-th entry in this vector is the function
%                     value associated with the i-th entry in the 'x'
%                     vector
%
%     output:
%             cspl    five column matrix containing the information which
%                     defines the "not-a-knot" cubic spline interpolant
%                     - first column: interpolating points
%                     - second column: function values
%                     - third column: coefficients of linear terms
%                     - fourth column: coefficients of quadratic terms
%                     - fifth column: coefficients of cubic terms
%
%     NOTE:
%             to evaluate the "not-a-knot" cubic spline interpolant apply
%             the routine SOLVE_SPLINE to the output of this routine 
%

n = length(x);
m = length(fx);
if n ~= m
    error('cubic_clamp:InvalidParameters', 'Length of parameters must match');
end

h = zeros(1, n - 1);
al = zeros(1, n);
l = al;
mu = l;
z = l;
b = l;
c = l;
d = l;

h(1) = x(2) - x(1);
al(1) = 3.0 * ( fx(2) - fx(1) ) / h(1) - 3.0 * fpO;
l(1) = 2 * h(1);
mu(1) = 0.5;
z(1) = al(1) / l(1);

for i = 2 : n - 1
    h(i) = x(i + 1) - x(i);
    al(i) = ( 3.0 / h(i) ) * ( fx(i + 1) - fx(i) ) - ( 3.0 / h(i - 1) ) * (fx(i) - fx(i - 1) );
    l(i) = 2.0 * ( x(i + 1) - x(i - 1) ) - h(i - 1) * mu(i - 1);
    mu(i) = h(i) / l(i);
    z(i) = ( al(i) - h(i - 1) * z(i - 1) ) / l(i);
end

al(n) = 3.0 * fpN - 3.0 * ( fx(n) - fx(n - 1) ) / h(n - 1);
l(n) = h(n - 1) * ( 2.0 - mu(n - 1) );
z(n) = ( al(n) - h(n - 1) * z(n - 1) ) / l(n);
b(n) = NaN;
c(n) = z(n);
d(n) = NaN;

for j = n - 1 : -1 : 1
    c(j) = z(j) - mu(j) * c(j + 1);
    b(j) = ( fx(j + 1) - fx(j) ) / h(j) - h(j) * ( c(j + 1) + 2.0 * c(j) ) / 3.0;
    d(j) = ( c(j + 1) - c(j) ) / ( 3.0 * h(j) );
end

cspl = [x' fx' b' c' d'];
end

function y = solve_spline( cspl, x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

end

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
