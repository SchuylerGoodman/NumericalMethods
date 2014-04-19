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

tic

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

toc
end

