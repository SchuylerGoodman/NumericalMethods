function cspd = diff_spline( cspl )
%DIFF_SPLINE Compute the derivative of the cubic spline interpolant
%
%     calling sequences:
%             cspd = diff_spline(cspl)
%
%     inputs:
%             cspl    five column matrix containing the information which
%                     defines the "not-a-knot" cubic spline interpolant
%                     - first column: interpolating points
%                     - second column: function values
%                     - third column: coefficients of linear terms
%                     - fourth column: coefficients of quadratic terms
%                     - fifth column: coefficients of cubic terms
%                     Typically this is the output of CUBIC_NAT,
%                     CUBIC_CLAMP, or DIFF_SPLINE
%
%     outputs:
%             cspd    five column matrix containing the information which
%                     defines the derivative of the "not-a-knot" cubic 
%                     spline interpolant
%                     - first column: interpolating points
%                     - second column: function values
%                     - third column: coefficients of linear terms
%                     - fourth column: coefficients of quadratic terms
%                     - fifth column: coefficients of cubic terms

%   Detailed explanation goes here
tic
[n m] = size(cspl);
cspd = zeros(n, m);
cspd(:, 1) = cspl(:, 1);
cspd(:, 2) = cspl(:, 3);
cspd(:, 3) = cspl(:, 4) * 2;
cspd(:, 4) = cspl(:, 5) * 3;
toc
end
