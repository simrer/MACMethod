function [f,grad] = AckleyFcn(x)
d = length(x); a = 20; b = 0.2; c = 2*pi; % recommended from litrature
sum1 = 0;
sum2 = 0;
sum3 = 0;
grad = zeros(1,d); 
for i = 1:d
    xi = x(i);
    sum1 = sum1 + xi^2;
    sum2 = sum2 + cos(c*xi);
    sum3 = sum3 + sin(c*xi);
end
term1 = -a*exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);
f = term1 + term2 + a + exp(1);
for i = 1:d
    xi = x(i);
    if nargout > 1
        
         grad(1,i) = term1*(-b*xi/(d*sqrt(sum1/d))) + ...
             term2*(-c/d*sum3); 
    end
end
end