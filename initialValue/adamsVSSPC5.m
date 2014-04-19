function [ t, yt ] = adamsVSSPC5( fty, lowerLimit, upperLimit, initialCondition, tolerance, stepmin, stepmax )
%ADAMSVSSPC5 Estimate value of the initial-value problem using the Adams
%               Variable Step-Size Predictor-Corrector Method. Uses the
%               Adams-Bashforth 5-step Method and the Adams-Moulton 4-step
%               Method for the Predictor-Corrector.
%
%           Inputs:
%               fty - f(t, y), or y'.  The differential equation whose
%                   solution is being estimated.
%               lowerLimit - Lower bound for t
%               upperLimit - Upper bound for t, value at which to estimate
%                   y(t)
%
%           Outputs:
%                   t - 
%

step = stepmax;
i = 1;

FLAG = 1;
LAST = 0;

disp('i     t              step         w');
fprintf('0     %6.5f       %6.3f        %6.8f\n', lowerLimit, step, initialCondition);
% use Runge-Kutta Order 4 to get next 4 values
[t, yt] = subrk4(fty, step, i, i + 4, lowerLimit, initialCondition);
NFLAG = 1;
% update i to account for Runge-Kutta values
i = i + 4; % 5

% we now estimate for next i
i = i + 1; % 6
% try this value of t
tt = t(4) + step;

while FLAG && abs(step) ~= 0
    % Adams-Bashforth 5-step Method to get an estimate
    wp = yt(i - 1) + step / 720 * ( 1901 * feval(fty, t(i - 1), yt(i - 1)) - 2774 * feval(fty, t(i - 2), yt(i - 2)) + 2616 * feval(fty, t(i - 3), yt(i - 3)) - 1274 * feval(fty, t(i - 4), yt(i - 4)) + 251 * feval(fty, t(i - 5), yt(i - 5)));
    % Adams-Moulton 3-step Method to correct the estimate
    wc = yt(i - 1) + step / 720 * ( 251 * feval(fty, tt, wp) + 646 * feval(fty, t(i - 1), yt(i - 1)) - 264 * feval(fty, t(i - 2), yt(i - 2)) + 106 * feval(fty, t(i - 3), yt(i - 3)) - 19 * feval(fty, t(i - 4), yt(i - 4)));
    % estimate error
    alpha = 19 * abs( wc - wp ) / ( 270 * step );
    
    if alpha <= tolerance % if the estimate is within the tolerance
        t = [ t tt ]; % accept this t
        yt = [ yt wc ]; % accept this y(t)
        
        if NFLAG
            for j = i - 4 : i
                fprintf('%d     %6.5f       %6.3f        %6.8f\n', j - 1, t(j), step, yt(j));
            end
        else
            fprintf('%d     %6.5f       %6.3f        %6.8f\n', i - 1, t(i), step, yt(i));
        end
        
        if LAST % if we calculated the last value
            FLAG = 0; % leave while loop and end computation
        else % otherwise
            i = i + 1; % calculating next i
            NFLAG = 0; % we accept all previous values
            
            % if estimate is too accurate or next step is outside bounds
            if alpha <= 0.1 * tolerance || t(i - 1) + step > upperLimit
                % get step size multiplier
                q = ( tolerance / ( 2 * alpha ) )^(1/4);
                % get new step size
                if q > 4 % upper bound is 4
                    step = 4 * step;
                else
                    step = q * step;
                end
                
                if step > stepmax % additional upper bound on stepmax
                    step = stepmax;
                end
                
                if t(i - 1) + 5 * step > upperLimit % if we're at the end
                    % make sure step size does not go out of bounds
                    step = ( upperLimit - t(i - 1) ) / 5;
                    LAST = 1;
                end
                % calculate next Runge-Kutta values
                [t, yt] = subrk4(fty, step, i - 1, i + 3, t, yt);
                NFLAG = 1; % we have not accepted these values yet
                i = i + 4; % go one past all calculated values
            end
        end
    else % if estimate is not within the tolerance
        % we need a new step size
        q = ( tolerance / ( 2 * alpha ) )^(1/4);
        if q < 0.1
            step = 0.1 * step;
        else
            step = q * step;
        end
        
        if step < stepmin
            FLAG = 0;
            disp('hmin exceeded.');
        else % if step size is OK
            if NFLAG % if we have not accepted the values
                i = i - 4; % recalculate last ones
            end
            % get next Runge-Kutta values
            [t, yt] = subrk4(fty, step, i - 1, i + 3, t, yt);
            i = i + 4;
            NFLAG = 1; % values not yet accepted
        end
    end
    tt = t(i - 1) + step; % next temporary t
end

end

function [ newt, newyt ] = subrk4(fty, step, initialI, finalI, t, yt)

assert(initialI <= length(t), 'Invalid parameters: initialI outside bounds.');
assert(finalI >= length(t), 'Cannot change values at given inputs.');

tLen = length(t);

newt = zeros(1, finalI);
newt(1 : tLen) = t(:);

newyt = newt;
newyt(1 : tLen) = yt(:);

for i = initialI + 1 : finalI
    newt(i) = newt(i - 1) + step;
    % yi = yi-1 + h / 4 * (f(ti-1, yi-1) + 3 * f(ti-1 + (2h / 3), wi-1 + (2h / 3) * f(ti-1 + h / 3, yi-1 + h/3 * f(ti-1, wi-1))))
    a = step * feval(fty, newt(i - 1), newyt(i - 1));
    b = step * feval(fty, newt(i - 1) + step / 2, newyt(i - 1) + 1 / 2 * a);
    c = step * feval(fty, newt(i - 1) + step / 2, newyt(i - 1) + 1 / 2 * b);
    d = step * feval(fty, newt(i), newyt(i - 1) + c);
    newyt(i) = newyt(i - 1) + 1 / 6 * (a + 2 * b + 2 * c + d);
end

end
