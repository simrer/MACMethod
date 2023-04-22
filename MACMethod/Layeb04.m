function f = Layeb04(x)
% multimodal, non-separable, scalable
% optimal value is (ln(0.001)-1)(dim-1) (-31.6310) and it occurs
% at alternation of 0 and (2k-1)Pi, k is an integer
% not reproduceble
% optimal value is f(x*) = (ln(0.001)-1)(dim-1) on the search space
% (-10,10)
n = length(x);
indx = zeros(1,n-1);
for i = 1:n-1
    indx(i) = log(abs(x(i)*x(i+1))+0.001)+cos(x(i)+x(i+1));
end
f = sum(indx);
