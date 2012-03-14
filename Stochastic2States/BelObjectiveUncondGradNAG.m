function [ grad, iflag] = BelObjectiveUncondGradNAG(n,x,iflag)
global V Vcoef R u2btild Par s_
%BELOBJECTIVEUNCOND Computes the Bellman objective with 
%   Detailed explanation goes here
    gamma = Par.gamma;
    sigma = Par.sigma;
    beta =  Par.beta;
    P = Par.P;
    theta_1 = Par.theta(1);
    theta_2 = Par.theta(2);
    g = Par.g;
    alpha = Par.alpha;
    if min(x)>0
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
    %compute components from unconstrained guess
    [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,sigma,gamma);
    [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,sigma,gamma,beta,P);
    
    %compute objective
    grad1 = zeros(3,1);
    X = [c2_1^(-sigma)*btildprime(1),c2_1^(-sigma)/c1_1^(-sigma)];%state next period
    Vobj = P(s_,1)*(alpha(1)*u(c1_1,l1(1),sigma,gamma)+alpha(2)*u(c2_1,l2(1),sigma,gamma)...
            +beta*funeval(Vcoef{1},V(1),X));
    dV = funeval(Vcoef{1},V(1),X,eye(2));
    grad1(1) = P(s_,1)*(alpha(1)*c1_1^(-sigma)+beta*sigma*c1_1^(sigma-1)*dV(2));
    grad1(2) = 0;
    grad1(3) = P(s_,1)*(alpha(2)*c2_1^(-sigma)-beta*sigma*c2_1^(-sigma-1)*(btildprime(1)+dV(2)));
    
    grad1 = grad1+P(s_,1)*( grad_btildprime(:,1)*c2_1^(-sigma)*beta*dV(1) - alpha(1)*l1(1)^gamma*l1grad(:,1)...
        -alpha(2)*l2(1)^gamma*l2grad(:,1));
    
    grad2 = zeros(3,1);
    X = [c2_2^(-sigma)*btildprime(2),c2_2^(-sigma)/c1_2^(-sigma)];%state next period
    Vobj = Vobj + P(s_,2)*(alpha(1)*u(c1_2,l1(2),sigma,gamma)+alpha(2)*u(c2_2,l2(2),sigma,gamma)...
            +beta*funeval(Vcoef{2},V(2),X));
    dV = funeval(Vcoef{2},V(2),X,eye(2));
    grad2(1) = 0;
    grad2(2) = P(s_,2)*(alpha(1)*c1_2^(-sigma)+beta*sigma*c1_2^(sigma-1)*dV(2));
    grad2(3) = 0;
    
    d_c2_2 = P(s_,2)*(alpha(2)*c2_2^(-sigma)-beta*sigma*c2_2^(-sigma-1)*(dV(2)+btildprime(2)));
    
    grad2 = grad2 + d_c2_2*grad_c2_2;
    grad2 = grad2 + P(s_,2)*( grad_btildprime(:,2)*c2_2^(-sigma)*beta*dV(1)-alpha(1)*l1(2)^gamma*l1grad(:,2)...
        -alpha(2)*l2(2)^gamma*l2grad(:,2));
    
    grad = grad1+grad2;
    if ~isreal(grad)
    
    grad=abs(grad)+100;
    end
    else
        grad=abs(x)+100;
    end

end


function [ c2_2 grad ] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma)

    %Compute c2_2 from formula
    frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) );
    c2_2 = frac^(-1/sigma);
    grad=zeros(3,1);
    %compute the gradients for c1_1,c1_2,c2_1
    grad(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2));
    grad(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R;
    grad(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2);
end

function [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,sigma,gamma)

    %Compute l1 form formula
    l1_1den = (theta_1+theta_2*(theta_2*c1_1^sigma/(theta_1*c2_1^sigma))^(1/gamma));
    l1(1) = (c1_1+c2_1+g(1))/l1_1den;
    l1_2den = (theta_1+theta_2*(theta_2*c1_2^sigma/(theta_1*c2_2^sigma))^(1/gamma));
    l1(2) = (c1_2+c2_2+g(2))/l1_2den;
    
    %compute gradients of l1(1) for c1_1,c1_2,c2_1
    l1grad(1,1) = (l1_1den-(c1_1+c2_1+g(1))*theta_2*sigma/gamma*...
        (theta_2*c1_1^sigma/(theta_1*c2_1^sigma))^(1/gamma)/c1_1)/l1_1den^2;
    l1grad(2,1) = 0;
    l1grad(3,1) = (l1_1den+(c1_1+c2_1+g(1))*theta_2*sigma/gamma*...
        (theta_2*c1_1^sigma/(theta_1*c2_1^sigma))^(1/gamma)/c2_1)/l1_1den^2;
    
    %compute gradients of l1(1) for c1_1,c1_2,c2_1
    l1grad(1,2) = 0;
    l1grad(2,2) = (l1_2den-(c1_2+c2_2+g(2))*theta_2*sigma/gamma*...
        (theta_2*c1_2^sigma/(theta_1*c2_2^sigma))^(1/gamma)/c1_2)/l1_2den^2;
    l1grad(3,2) = 0;
    %use chain rule for c2_2
    d_c2_2 = (l1_2den+(c1_2+c2_2+g(2))*theta_2*sigma/gamma*...
        (theta_2*c1_2^sigma/(theta_1*c2_2^sigma))^(1/gamma)/c2_2)/l1_2den^2;
    l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2;
    
    %compute l2 from formula
    l2_1den = (theta_1*(theta_1*c2_1^sigma/(theta_2*c1_1^sigma))^(1/gamma)+theta_2);
    l2(1) = (c1_1+c2_1+g(1))/l2_1den;
    l2_2den = (theta_1*(theta_1*c2_2^sigma/(theta_2*c1_2^sigma))^(1/gamma)+theta_2);
    l2(2) = (c1_2+c2_2+g(2))/l2_2den;
    
    %compute gradients of l2(1) for c1_1,c1_2,c2_1
    l2grad(1,1) = ( l2_1den+(c1_1+c2_1+g(1))*theta_1*sigma/gamma*...
        (theta_1*c2_1^sigma/(theta_2*c1_1^sigma))^(1/gamma)/c1_1)/l2_1den^2;
    l2grad(2,1) = 0;
    l2grad(3,1) = ( l2_1den-(c1_1+c2_1+g(1))*theta_1*sigma/gamma*...
        (theta_1*c2_1^sigma/(theta_2*c1_1^sigma))^(1/gamma)/c2_1)/l2_1den^2;
    
    %compute gradients of l2(2) for c1_1,c1_2,c2_1
    l2grad(1,2) = 0;
    l2grad(2,2) = (l2_2den+(c1_2+c2_2+g(2))*theta_1*sigma/gamma*...
        (theta_1*c2_2^sigma/(theta_2*c1_2^sigma))^(1/gamma)/c1_2)/l2_2den^2;
    l2grad(3,2) = 0;
    %use chain rule to get the effect of c2_2
    d_c2_2 = (l2_2den-(c1_2+c2_2+g(2))*theta_1*sigma/gamma*...
        (theta_1*c2_2^sigma/(theta_2*c1_2^sigma))^(1/gamma)/c2_2)/l2_2den^2;
    l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
    
end

function [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,sigma,gamma,beta,P)
    %get expected value of marginal utility of agent 2
    Eu2 = P(s_,1)*c2_1^(-sigma)+P(s_,2)*c2_2^(-sigma);
    
    %compute btildeprime from formula
    btildprime(1) = u2btild/(beta*Eu2)...
        +c1_1-c2_1-l1(1)^(1+gamma)/c1_1^(-sigma)+l2(1)^(1+gamma)/c2_1^(-sigma);

    %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
    grad_btildprime(1,1) = 1-sigma*l1(1)^(1+gamma)*c1_1^(sigma-1);
    grad_btildprime(2,1) = 0;
    grad_btildprime(3,1) =sigma*u2btild*P(s_,1)*c2_1^(-sigma-1)/(beta*Eu2^2)...
        -1+sigma*l2(1)^(1+gamma)*c2_1^(sigma-1);

    %figure out their affects through c2_2, l1_1,l2_1
    d_c2_2 = sigma*u2btild*P(s_,2)*c2_2^(-sigma-1)/(beta*Eu2^2);
    d_l1_1 = -(1+gamma)*l1(1)^gamma*c1_1^sigma;
    d_l2_1 = (1+gamma)*l2(1)^gamma*c2_1^sigma;
    grad_btildprime(:,1) = grad_btildprime(:,1) + d_c2_2*grad_c2_2+d_l1_1*l1grad(:,1)+d_l2_1*l2grad(:,1);

    %Compute btildprime(2) from formula
    btildprime(2) = u2btild/(beta*Eu2)...
        +c1_2-c2_2-l1(2)^(1+gamma)/c1_2^(-sigma)+l2(2)^(1+gamma)/c2_2^(-sigma);

    %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
    grad_btildprime(1,2) = 0;
    grad_btildprime(2,2) = 1-sigma*l1(2)^(1+gamma)*c1_2^(sigma-1);
    grad_btildprime(3,2) = sigma*u2btild*P(s_,1)*c2_1^(-sigma-1)/(beta*Eu2^2);
    %figure out their affects through c2_2, l1_2,l2_2
    d_c2_2 = sigma*u2btild*P(s_,2)*c2_2^(-sigma-1)/(beta*Eu2^2)-1+sigma*l2(2)^(1+gamma)*c2_2^(sigma-1);
    d_l1_2 = -(1+gamma)*l1(2)^gamma*c1_2^sigma;
    d_l2_2 = (1+gamma)*l2(2)^gamma*c2_2^sigma;

    grad_btildprime(:,2) = grad_btildprime(:,2) + d_c2_2*grad_c2_2+d_l1_2*l1grad(:,2)+d_l2_2*l2grad(:,2);

end