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

tic
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
toc
end
