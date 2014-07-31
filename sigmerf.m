function s = sigmerf(x,G,P)

%s = 1+erf(P.*(x-G));
s = 1+tanh(P.*(x-G));

end