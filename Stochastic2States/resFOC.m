function  [ res, iflag] =resFOC(n,x,iflag)
global V Vcoef R u2btild Par s_ flagCons

theta_1=Par.theta_1;
sigma=Par.sigma;
gamma=Par.gamma;
g=Par.g;
alpha_1=Par.alpha_1;
alpha_2=Par.alpha_2;
P=Par.P;
beta=Par.beta;
u2btildLL=Par.u2btildLL;
u2btildUL=Par.u2btildUL;



%% Get the variables from x

c1(1)=x(1); %consumption of agent 1 state 1
c1(2)=x(2); %consumption of agent 1 state 2

c2(1)=x(3);  %consumption of agent 2 state 1
c2(2)=x(4); %consumption of agent 2 state 2

l1(1)=x(5); %labor supply of agent 1 state 1
l1(2)=x(6); %labor supply of agent 1 state 2

% Now get the values of u2btildprime. The following section uses the
% flagCons to see which of the bounds are binding. It accordingly set the
% value of u2btildprime to the limit and solves for the mutltiplier
MuU=zeros(1,2);
MuL=zeros(1,2);
switch flagCons
    case 'LL_'
       % lower limit binds for state 1 only
       MuL(1)=x(7);
       MuL(2)=0;
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=x(8);
       
    case '_LL'
       % lower limit binds for state 2 only
       MuL(1)=0;
       MuL(2)=x(8);
       u2btildprime(1)=x(7);
       u2btildprime(2)=u2btildLL;
       
    case 'LLLL'
      % lower limit binds for both the states
       MuL(1)=x(7);
       MuL(2)=x(8);
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=u2btildLL;     
        
        
        
    case 'UL_'
     % upper limit binds for state 1 only

       MuU(1)=x(7);
       MuU(2)=0;
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=x(8);
       
        
    case '_UL'
         % upper limit binds for state 2 only
       MuU(1)=0;
       MuU(2)=x(8);
       u2btildprime(1)=x(7);
       u2btildprime(2)=u2btildUL;
        
        
    case 'ULUL'
        
        
       % upper limit binds for both the states
       MuL(1)=x(7);
       MuL(2)=x(8);
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=u2btildUL;     
        
        
        
end


lambda_I(1)=x(9);  % Multiplier on I(1)
lambda_I(2)=x(10); % Multiplier on I(2)

lambda_B=x(11);  % Multiplier on B

lambda_R(1)=x(12); % Multiplier on R(1)
lambda_R(2)=x(13); % Multiplier on R(2)


% Check for non negativity of consumption and labor
if min(x(1:6))>0
    
    
    % Get the expected value of the marginal utilities 
    Eu1=P(s_,1)*c1(1)^(-sigma)+P(s_,2)*c1(2)^(-sigma);
        Eu2=P(s_,1)*c2(1)^(-sigma)+P(s_,2)*c2(2)^(-sigma);

% Derivative of the Implementability constraint
    dI(:,1)=[-1+l1(1)^(1+gamma)*sigma*c1(1)^(sigma-1); ...                                                               % c1(1)
         0; ...                                                                                                          % c1(2)
        1+u2btildprime(1)*sigma*c2(1)^(sigma-1)-((sigma*u2btild*P(s_,1)*c2(1)^(-sigma-1))/(beta*Eu2^2)); ....             %c2(1)
        -(u2btild*P(s_,2)*sigma*c2(2)^(-sigma-1)/(beta*Eu2^2));...                                                         %c2(2)
         (1+gamma)*l1(1)^gamma*c1(1)^(sigma);...                                                                         %l1(1)
         0; ...                                                                                                           %l1(2)
         c2(1)^(sigma); ...                                                                                             %u2btildprime(1)
         0;];                                                                                                           %u2btldprime (2)   
        
     
    dI(:,2)=[0; ...                                                                                                     %c1(1)
        -1+l1(2)^(1+gamma)*sigma*c1(2)^(sigma-1); ...                                                                   %c1(2)
        -(u2btild*P(s_,1)*sigma*c2(1)^(-sigma-1)/(beta*Eu2^2));   ...                                                     %c2(1)
        1+u2btildprime(2)*sigma*c2(2)^(sigma-1)-((sigma*u2btild*P(s_,2)*c2(2)^(-sigma-1))/(beta*Eu2^2)); ....            %c2(2)
        0;                                                                                                              %l1(1)
       (1+gamma)*l1(2)^gamma*c1(2)^(sigma);...                                                                          %l1(2)
          0;                                                                                                            %u2btildprime(1)
         c2(2)^(sigma);];                                                                                            %u2btildprime(2)
         
     
     % Derivative of the Bond Pricing constraint
 dB=[ (Eu2/Eu1^2)*P(s_,1)*sigma*c1(1)^(-sigma-1);...                                            %c1(1)
      (Eu2/Eu1^2)*P(s_,2)*sigma*c1(2)^(-sigma-1);...                                            %c1(2)
      (-1/Eu1)*P(s_,1)*sigma*c2(1)^(-sigma-1);...                                               %c2(1)
      (-1/Eu1)*P(s_,2)*sigma*c2(2)^(-sigma-1);...                                               %c2(2)
      0; ...                                                                                  %l1(1)
     0; ...                                                                                   %l1(2)
     0; ...                                                                                   %u2btildprime(1)
     0; ...                                                                                   %u2btildprime(2)
     ];
 
 % Derivative of the resource constraint
 dR=[1 0; ...                                                                                    %c1(1)
     0 1; ...                                                                                    %c1(2)
     1 0; ...                                                                                    %c2(1)
     0 1; ...                                                                                    %c2(2)
     -theta_1 0; ...                                                                             %l1(1)    
     0 -theta_1; ...                                                                             %l1(2)
     0 0; ...                                                                                    %u2btildprime(1)    
     0 0;];                                                                                      %u2btildprime(1) 
 
 % Derivative of the inequality constraints
 dLL=[0 0; ...                                                                                   %c1(1)
     0 0; ...                                                                                    %c1(2)
     0 0; ...                                                                                    %c2(1)
     0 0; ...                                                                                    %c2(2)
     0 0; ...                                                                                    %l1(1)    
     0 0; ...                                                                                    %l1(2)
     1 0; ...                                                                                    %u2btildprime(1)    
     0 1;];                                                                                      %u2btildprime(1) 
 
 dUL=[0 0; ...                                                                                   %c1(1)
     0 0; ...                                                                                    %c1(2)
     0 0; ...                                                                                    %c2(1)
     0 0; ...                                                                                    %c2(2)
     0 0; ...                                                                                    %l1(1)    
     0 0; ...                                                                                    %l1(2)
     -1 0; ...                                                                                   %u2btildprime(1)    
     0 -1;];                                                                                     %u2btildprime(1) 
 
%% Gradient of the value function 

% tomorrow
X(1,:)=[u2btildprime(1) c2(1)^(-sigma)/c1(1)^(-sigma)];
X(2,:)=[u2btildprime(2) c2(2)^(-sigma)/c1(2)^(-sigma)];
dV_u2b(1)=funeval(Vcoef{1},V(1),X(1,:),[1 0 ]);
dV_u2b(2)=funeval(Vcoef{2},V(2),X(2,:),[ 1 0]);
dV_R(1)=funeval(Vcoef{1},V(1),X(1,:),[0 1 ]);
dV_R(2)=funeval(Vcoef{2},V(2),X(2,:),[0 1 ]);


% Derivative of the objective constraint
dVobj(:,1) =[ alpha_1*P(s_,1)*(c1(1)^(-sigma)+beta*dV_R(1)*sigma*c1(1)^(sigma-1)*c2(1)^(-sigma));...        %c1(1)
            0; ...                                                                                       %c1(2)     
            alpha_2*P(s_,1)*(c2(1)^(-sigma)+beta*dV_R(1)*c1(1)^(sigma)*c2(1)^(-sigma-1)*(-sigma)); ...   %c2(1)
            0;                                                                                           %c2(2)    
            -alpha_1*P(s_,1)*l1(1)^(gamma); ...                                                          %l1(1)
            0;  ...                                                                                         %l1(2)    
            beta*P(s_,1)*dV_u2b(1); ....                                                                 %u2btildeprime(1)    
            0;                      ...                                                                  %u2btildprime(2)
            ];
        

dVobj(:,2) =[ 0; ...                                                                                       %c1(1)    
            alpha_1*P(s_,2)*(c1(2)^(-sigma)+beta*dV_R(2)*sigma*c1(2)^(sigma-1)*c2(2)^(-sigma));...        %c1(2)    
            0; ...                                                                                       %c2(1)    
            alpha_2*P(s_,2)*(c2(2)^(-sigma)+beta*dV_R(2)*c1(2)^(sigma)*c2(2)^(-sigma-1)*(-sigma)); ...   %c2(2)
            0; ...                                                                                       %l1(1)
            -alpha_1*P(s_,2)*l1(2)^(gamma); ...                                                          %l1(2)
            0;                                                                                           %u2btildprime(1)    
            beta*P(s_,2)*dV_u2b(2); ....                                                                 %u2btildprime(2)
            ];                                                                                          
        
        
        
FOC=dVobj(:,1)+dVobj(:,2) +lambda_I(:,1)*dI(:,1) + lambda_I(:,2)*dI(:,2) + lambda_B*dB +lambda_R(:,1)*dR(:,1) + lambda_R(:,2)*dR(:,2) + MuL(:,1)*dLL(:,1)+MuL(:,2)*dLL(:,2)+ MuU(:,1)*dUL(:,1)+MuU(:,2)*dUL(:,2);            


EqCons=[(c2(1)-c1(1))+l1(1)^(1+gamma)/c1(1)^(-sigma)+u2btildprime(1)/c2(1)^(-sigma)-u2btild/(beta*Eu2);...         % Implementability state 1
    (c2(2)-c1(2))+l1(2)^(1+gamma)/c1(2)^(-sigma)+u2btildprime(2)/c2(2)^(-sigma)-u2btild/(beta*Eu2);...             % Implementability state 2
Eu2/Eu1-R; ...                                                                                   % Bond Pricing   
c1(1)+c2(1)+g(1)-theta_1*l1(1); ...                                                               % Resource Constraint state 1
c1(2)+c2(2)+g(2)-theta_1*l1(2);]; 


res=[FOC ;EqCons];
else
    %disp(x)
    res=abs(x)+100;

end
