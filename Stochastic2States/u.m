function u=u(c,l,sigma,gamma)

if sigma~=1 
uc=c.^(1-sigma)/(1-sigma);
else
    uc=log(c);
end
if gamma~=-1
ul=-l.^(1+gamma)/(1+gamma);
else
    ul=log(l);
end
u=uc+ul;