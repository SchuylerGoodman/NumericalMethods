function [p]=Newton1(f,df,p0,tol)
tic;
while 1
    %df=(f(p0+tol)-f(p0))/tol;
p=p0-f(p0)/df(p0);
if abs(p-p0)<tol, break; end
p0=p;
end
toc