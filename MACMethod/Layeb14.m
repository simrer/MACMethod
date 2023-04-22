function y = Layeb14(x)
% multimodal, non-separable and scalable function. It is very hard
% to optimize and most algorithms fail. search domain is (-100,100)
% optimal value is 0 and occurs at permutation of 0 and -1
n = length(x);
indx = zeros(1,n-1);
for i = 1:n-1
    indx(i) = 100*abs(x(i)^2-x(i+1))^0.1+abs(log10((x(i)+x(i+1)^2)));
end
y = -sum(indx);