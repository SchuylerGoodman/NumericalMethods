classdef rootfinder < handle
    %ROOTFINDER - finds all roots for a given polynomial or function
    %   Accepts for input either row vector polynomials, function handles,
    %   or symbolic functions.
    
    properties
        
        tolerance                   % tolerance for estimation
        
        max_iter                    % maximum number of iterations
        
        current_f                   % current polynomial function (in row vector format)
        
        p0                          % x which is currently estimated
        
    end
    
    properties (SetAccess = 'private')
        
        X                           % row vector containing all roots
        
    end
    
    properties (GetAccess = 'private')
        
        f_less                      % divided polynomial function (in row vector format)
        
    end
    
    methods (Access = 'private')
        
        % find all roots
        function P = roots(self)
            self.f_less = self.current_f;
            self.f_less = self.checkPoly();
            len = length(self.current_f);
            P = zeros(1, len - 1);
            iter = 1;
            while len > 1
                P(iter) = self.root();
                if isnan(P(iter))
                    error('rootfinder:NaN', 'Error finding root.');
                end
                [self.f_less ~] = self.synthDivide(P(iter));
                len = length(self.f_less);
                iter = iter + 1;
            end
        end
        
        % find single root
        function p = root(self)
            p = self.newton();
            if isnan(p)
                p = self.bisection();
                if isnan(p)
                    p = self.mullers();
                end
            end
        end
        
        % newton's method
        function p = newton(self)
            disp('Trying Newton''s Method...');
            disp('i       x               err');
            iter = 0;
            p = NaN;
            p_naught = self.p0;
            
            % calculate lower max iterations, so we can switch to modified
            % if necessary (I will have to adjust this, I think...)
            max_iter_newton = 10; % self.max_iter / (self.max_iter / 10);
            while iter < max_iter_newton
                % estimate df
                df = (polyval(self.f_less, p_naught + self.tolerance) - ...
                    polyval(self.f_less, p_naught)) / self.tolerance;
            
                % if df is close enough to equal zero, try another method
                if abs(df) < self.tolerance
                    p = NaN;
                    return;
                end
                
                % iterate p
                p = p_naught - polyval(self.f_less, p_naught) / df;
                fprintf('%d     %6.5f           %6.5f\n', iter, p, abs(p - p_naught));
                
                % we found it!
                if abs(p - p_naught) < self.tolerance
                    return;
                end
                iter = iter + 1;
                p_naught = p;
            end
            if iter == max_iter_newton
                p = self.newton_mod(p_naught);
            end
        end
        
        % modified newton's method
        function p = newton_mod(self, p0)
            disp('Possible not a simple zero - trying Modified Newton''s Method');
            disp('i       x               err');
            iter = 0;
            p_naught = p0;
            
            while iter < self.max_iter
                % estimate df
                df = (polyval(self.f_less, p_naught + self.tolerance) - ...
                    polyval(self.f_less, p_naught)) / self.tolerance;
                df2 = ( (df + self.tolerance) - (df) ) / self.tolerance;
            
                % if both df and df2 is close enough to equal zero, try another method
                if abs(df) < self.tolerance && abs(df2) < self.tolerance
                    p = NaN;
                    return;
                end
                
                % calculate mu and mu prime
                m = polyval(self.f_less, p_naught) / df;
                mp = (df^2 - polyval(self.f_less, p_naught) * df2) / df^2;
                
                % iterate p
                p = p_naught - m / mp;
                fprintf('%d     %6.5f           %6.5f\n', iter, p, abs(p - p_naught));
                
                % we found it!
                if abs(p - p_naught) < self.tolerance
                    return;
                end
                iter = iter + 1;
                p_naught = p;
            end
            p = NaN;
        end
        
        % bisection
        function p = bisection(self)
            disp('Trying Bisection method.');
            p = NaN;
            n = polyval(self.f_less, self.p0);
            iter = 0;
            
            getNegative = 0;
            if n < 0
                a = n;
            else
                b = n;
                getNegative = 1;
            end
            opp = self.getOpposite(getNegative);
            if isnan(opp)
                fprintf('In Bisection: Could not find f(x) with sign opposite to %6.5f.\n', n);
                p = opp;
                return;
            elseif getNegative
                a = opp;
            else
                b = opp;
            end
            
            % Now we have bounds for p
            disp('i       x               err');
            while 1  % convergence is guaranteed at this point (for continuous functions),
                % so loop until we converge.
                p = (a + b) / 2;
                if abs(p - a) < self.tolerance
                    return;
                elseif abs(p - a) > 1 / self.tolerance
                    p = NaN;  % it is probably not continuous
                    return;
                end
                f_a = polyval(self.f_less, a);
                %f_b = polyval(self.f_less, b);
                f_p = polyval(self.f_less, p);
                fprintf('%d     %6.5f           %6.5f\n', iter, p, abs(p - a));
                if f_a * f_p > 0
                    a = p;
                else
                    b = p;
                end
                iter = iter + 1;
            end
        end
        
        % method to get value with opposite sign for bisection
        % Uses various guesses to try to find a value.
        %
        % param negative - boolean (0 or 1) stating if the function needs
        % to return negative or positive (negative if 1, positive if 0).
        function oppositeValue = getOpposite(self, negative)
            tryInput = self.tolerance;
            oppositeValue = NaN;
            while tryInput <= 1 / self.tolerance
                tryValue = polyval(self.f_less, tryInput);
                tryInput = tryInput * 10;
                if (abs(tryValue) < self.tolerance) || (negative == (tryValue < 0))
                    oppositeValue = tryValue;
                    return;
                end
            end
        end
        
        % mullers method for complex zeros
        function p = mullers(self)
            disp('Complex zeros may exist - Trying Muller''s Method');
            disp('i       x                        err');
            p = NaN;
            iter = 3;
            p_naught = self.p0;
            p1 = p_naught;
            p2 = p_naught;
            multiplier = 1;

            while iter < self.max_iter
                while polyval(self.f_less, p_naught) == polyval(self.f_less, p1) && ...
                      polyval(self.f_less, p_naught) == polyval(self.f_less, p1)
                    multiplier = multiplier * 5;
                    if sign(p1) == -1
                        p_naught = p1 * multiplier;
                        p2 = p1 * multiplier * -1;
                    elseif sign(p1) == 1
                        p_naught = p1 * multiplier * -1;
                        p2 = p1 * multiplier;
                    else
                        p_naught = multiplier * -1;
                        p2 = multiplier;
                    end
                    if multiplier > 1 / self.tolerance
                        disp('This is probably a horizontal line.');
                        return;
                    end
                end
                
                h1 = p1 - p_naught;
                h2 = p2 - p1;
                d1 = (polyval(self.f_less, p1) - polyval(self.f_less, p_naught)) / h1;
                d2 = (polyval(self.f_less, p2) - polyval(self.f_less, p1)) / h2;
                d = (d2 - d1) / (h2 + h1);
                b = d2 + h2 * d;
                D = sqrt(b^2 - 4 * polyval(self.f_less, p2) * d);

                if abs(b - D) < abs(b + D)
                    E = b + D;
                else
                    E = b - D;
                end
                h = -2 * polyval(self.f_less, p2) / E;
                p = p2 + h;
                fprintf('%1d       %s         %6.5f\n', iter, num2str(p), abs(h));
                if abs(h) < self.tolerance
                    break;
                end
                p_naught = p1;
                p1 = p2;
                p2 = p;
                iter = iter + 1;
            end
        end
        
        % process function and return polynomial row vector
        function f = checkPoly(self)
            try
                polyval(self.current_f, self.p0);
            catch err
                if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                    error('MATLAB:rootfinder:invalid_input', ...
                        'Input function must be a row vector polynomial.');
                end
            end
            f = self.current_f;
        end
        
        % synthetic division for polynomial
        function [f remainder] = synthDivide(self, root)
            len = length(self.f_less);
            g = zeros(1, len - 1);
            next = 0;

            try
                polyval(self.f_less, 0);
            catch err
                error('MATLAB:syndiv:incorrectInput', ...
                    'Input must be a polynomial (row vector).');
            end

            if len > 1
                for i = 1 : len - 1
                    next = self.f_less(i) + ( root * next );
                    g(i) = next;
                end
            end
            remainder = self.f_less(len) + ( root * next );
            f = g;
        end
        
    end
    
    %======================================================================
    % PUBLIC METHODS
    %======================================================================
    methods
        
        % Constructor and Main
        function self = rootfinder(f, p0, tol, max_N)
            tic
            self.tolerance = tol;
            self.current_f = f;
            self.max_iter = max_N;
            self.p0 = p0;
            self.X = roots(self);
            disp(to_string(self));
            toc
        end
        
        % Output last function, p0, and estimated zeros
        function str = to_string(self)
            str = sprintf('Evaluation of %s at %d:\n%s', sprintf('%6.5f ', self.current_f), self.p0, sprintf('X = %s\n', num2str(self.X)));
        end
        
        % Solve with a different initial guess
%         function allX = solve(self, p0)
%             self.p0 = p0;
%             self.X = roots(self);
%             allX = self.X;
%         end
%         
        function allX = solve(self, f, p0)
            tic
            self.p0 = p0;
            self.current_f = f;
            self.X = roots(self);
            allX = self.X;
            toc
        end
    end
    
end

