function [xn yn] = parametric_herm(x, y, gx, gy)

%BEZIER     Constructs cubic Hermite polynomials with the given parameters
%           in parametric form.
%
%     inputs:
%             x       A ROW vector containing the x coordinates of all
%                     endpoints of the polynomials to construct.
%             y       A ROW vector containing the y coordinates of all
%                     endpoints of the polynomials to construct.
%             gx      A matrix containing the x coordinates of all
%                     guidepoints of the curves to construct.  The first
%                     column contains left guidepoints, and the second
%                     column contains right guidepoints.
%             gy      A matrix containing the y coordinates of all
%                     guidepoints of the curves to construct.  The first
%                     column contains left guidepoints, and the second
%                     column contains right guidepoints.
%
%     outputs:
%             xn      A matrix of coefficients of the constructed
%                     parametric polynomials that calculate xn(t). Each row
%                     represents x(t) for n = row.
%             yn      A matrix of coefficients of the constructed
%                     parametric polynomials that calculate yn(t). Each row
%                     represents y(t) for n = row.
%
%     notes:
%             IMPORTANT:  Input parameters must all be of the same length.

x_len = length(x);
y_len = length(y);
[gx_len ~] = size(gx);
[gy_len ~] = size(gy);

if x_len < 2 || x_len ~= y_len
    error('PARAMETRIC_HERM:InvalidEndpointInput', ...
        'First two parameters must have same length > 1');
end
if gx_len ~= x_len - 1 || gx_len ~= gy_len
    error('PARAMETRIC_HERM:InvalidGuidepointInput', ...
        'Last two parameters must have same length equal to first two parameters minus 1.');
end

xn = zeros(x_len - 1, 4);
yn = zeros(y_len - 1, 4);
left = 1;
right = 2;

for i = 1 : x_len - 1
    xn(i, 1) = x(i);
    xn(i, 2) = gx(i, left) - x(i);
    xn(i, 3) = 2 * x(i + 1) - x(i) - 2 * gx(i, left) + gx(i, right);
    xn(i, 4) = x(i) - x(i + 1) - gx(i, right) + gx(i, left);
    
    yn(i, 1) = y(i);
    yn(i, 2) = gy(i, left) - y(i);
    yn(i, 3) = 2 * y(i + 1) - y(i) - 2 * gy(i, left) + gy(i, right);
    yn(i, 4) = y(i) - y(i + 1) - gy(i, right) + gy(i, left);
end 

end
