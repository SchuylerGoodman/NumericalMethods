function [] = plot_param(xn, yn, t_low, t_high)

%PLOT_PARAM  Plots a series of parametric functions from t_low 
%            to t_high
%
%      inputs:
%              xn       A matrix containing x(t) CUBIC parametric 
%                       functions. Each row contains the coefficients in
%                       reverse order. Number of rows MUST BE the same as
%                       the number of rows in yn.
%                       (i.e. coeff1 + coeff2*t + coeff3*t^2 + coeff4*t^3)
%              yn       A matrix containing y(t) CUBIC parametric 
%                       functions. Each row contains the coefficients in
%                       reverse order. Number of rows MUST BE the same as
%                       the number of rows in xn.
%                       (i.e. coeff1 + coeff2*t + coeff3*t^2 + coeff4*t^3)
%              t_low    Lower bound for t for plotting functions.
%              t_high   Upper bound for t for plotting functions.
%
%      outputs:
%              All functions in xn and yn plotted from t_low to t_high.
%
%      NOTE:
%              Typically used with output of PARAMETRIC_HERM or BEZIER
%

[x_len ~] = size(xn);
[y_len ~] = size(yn);

if x_len < 1 || x_len ~= y_len
    error('PLOT_PARAM:InvalidParameters', ...
        'Input parameters must be same length.');
end
if t_low >= t_high
    error('PLOT_PARAM:InvalidTRange', ...
        't_low must be strictly less than t_high');
end

syms t;
ti = linspace(t_low, t_high, 1000);
hold on;

for i = 1 : x_len
    xt = 0;
    yt = 0;
    xi_len = length(xn(i, :));
    yi_len = length(yn(i, :));
    
    if xi_len ~= 4 || xi_len ~= yi_len
        error('PLOT_PARAM:InvalidParameters', ...
            'Input parameters must contain equal number of row vectors of length 4.');
    end
    
    for j = 1 : xi_len
        xt = xt + xn(i, j).* ti.^(j - 1);
        yt = yt + yn(i, j).* ti.^(j - 1);
    end
    
    plot(xt, yt);
    
end
hold off;

end

