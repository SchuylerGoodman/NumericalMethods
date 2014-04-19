function [app] = adaptiveQuad(f, left, right, tol, maxLevel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

level = 0;

step = ( right - left) / 2;
fLeft = feval(f, left);
fMid = feval(f, left + step);
fRight = feval(f, right);

% Make first approximation
parentValue = ( step / 3 ) * ( fLeft + 4 * fMid + fRight);

% Start recursion
app = subdivide(f, parentValue, left, left + step, right, step / 2, 10 * tol, level + 1, maxLevel);

end

function [sub] = subdivide(f, parentValue, left, middle, right, step, tol, level, maxLevel)

if level > maxLevel
    error('Output level %d exceeded', maxLevel);
end

fLeft = feval(f, left);
fLeftMid = feval(f, left + step);
fMid = feval(f, middle);
fRightMid = feval(f, middle + step);
fRight = feval(f, right);

% Compute subdivisions
leftValue = ( step / 3 ) * (fLeft + 4 * fLeftMid + fMid);
rightValue = ( step / 3 ) * (fMid + 4 * fRightMid + fRight);

% If our subdivisions agree with the parent to within the tolerance
if abs(leftValue + rightValue - parentValue) < tol
    sub = leftValue + rightValue;
    return; % Return the sum
else % If they do not agree, continue to subdivide
    leftValue = subdivide(f, leftValue, left, left + step, middle, step / 2, tol / 2, level + 1, maxLevel);
    rightValue = subdivide(f, rightValue, middle, middle + step, right, step / 2, tol / 2, level + 1, maxLevel);
end

% Return the sum
sub = leftValue + rightValue;
return;

end