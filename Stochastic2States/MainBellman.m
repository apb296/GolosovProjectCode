% This is the main file for executing the Bellman Equation
clc
close all
clear all
%%
% This script sets up the para structure and records a tex table with the
% parameters
SetParaStruc

%% Build Grid for the state variables
% This setups up the functional space and the grid. It also creates a table
% with the information about the grid and 2 3-d pictures visually showing
% the grid space for each value of the discrete state
BuildGrid

%% Set the Parallel Config
err='';
try
    matlabpool
catch err
end
if ~strcmpi(err.identifier,'MATLAB:UndefinedFunction')

if(matlabpool('size') > 0)
   matlabpool close
end

matlabpool open local;
end

%% Computing the  Initial Guess.
% This section uses the static first best without taxes to compute the
% initial guess for the value function coeffecients
for s_=1:sSize
    n=1;
    c1_=.1;
    for Rctr=1:RGridSize
        
        for u2btildctr=1:u2btildGridSize
            
            g_=g(s_);
            u2btild_=u2btildGrid(u2btildctr);
            R_=RGrid(Rctr);
            x_state_(s_,n,:)=[u2btild_ R_ ];
            R_0=(R_)^(-sigma);
            u1btild=u2btild_/R_;
            
            [c1_,fval,exitflagB,output]= fzero(@(c1) (R_0-1)*c1^(1-sigma)+((c1*(1+R_0)+g_)/theta_1)^(1+gamma)+u1btild*(1-1/beta),c1_);
            c1(s_)=c1_;
            c2(s_)=R_0*c1_;
            u2btildPrime(s_)=u1btild*R_;
            l1(s_)=(c1_+R_0*c1_+g_)/theta_1;
            % DetData=load('C:\Users\anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Deterministic1State\Data\c1000.mat');
            % cd ('C:\Users\anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Deterministic1State\')
            %  [PolicyRulesInit]=GetInitialApproxPolicy(u1btild,DetData.x_state,DetData.PolicyRulesStore);
            % [PolicyRules, V_new,exitflag,exitflagB]=SolveGradNAG(u1btild,DetData.c',DetData.V,PolicyRulesInit(1),DetData.Para) ;
            % cd ('C:\Users\anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Stochastic2States')
            % VDet(s_,n)=funeval(DetData.c,DetData.V,u1btild);
            
            V0(s_,n)=(alpha_1*u(c1(s_),l1(s_),sigma,gamma)+alpha_2*u(c2(s_),0,sigma,gamma))/(1-beta);
            xInit_0(s_,n,:)=[c1(s_) c2(s_) l1(s_) u2btildPrime(s_)/c2(s_)^(-sigma) u2btildPrime(s_) R_];
            n=n+1;
        end
    end
    c0(s_,:)=funfitxy(V(s_),squeeze(x_state_(s_,:,:)),squeeze(V0(s_,:))' );
    
    
end
x_state=vertcat([squeeze(x_state_(1,:,:)) ones(length(x_state_),1)] ,[squeeze(x_state_(1,:,:)) 2*ones(length(x_state_),1)]);
c=c0;

% slicing the state space for parfor loop later
u2btild_slice=x_state(:,1) ;
R_slice=x_state(:,2) ;
s_slice=x_state(:,3) ;
% This stores the values of the policy functions and multipliers that last
% worked

PolicyRulesWorked=[xInit_0(1,1,1) xInit_0(2,1,1) xInit_0(1,1,2)];

% This stores the policy rules for each point in the state
% space.
PolicyRulesStore1=[squeeze(xInit_0(1,:,1))' squeeze(xInit_0(1,:,1))' ...
    squeeze(xInit_0(1,:,2))' squeeze(xInit_0(1,:,2))'...
    squeeze(xInit_0(1,:,3))' squeeze(xInit_0(1,:,3))' ...
    squeeze(xInit_0(1,:,4))' squeeze(xInit_0(1,:,4))' ....
    squeeze(xInit_0(1,:,5))' squeeze(xInit_0(1,:,5))' ....
    squeeze(xInit_0(1,:,6))' squeeze(xInit_0(1,:,6))' ....
    ];
PolicyRulesStore2=[squeeze(xInit_0(2,:,1))' squeeze(xInit_0(2,:,1))' ...
    squeeze(xInit_0(2,:,2))' squeeze(xInit_0(2,:,2))'...
    squeeze(xInit_0(2,:,3))' squeeze(xInit_0(2,:,3))' ...
    squeeze(xInit_0(2,:,4))' squeeze(xInit_0(2,:,4))' ....
    squeeze(xInit_0(2,:,5))' squeeze(xInit_0(2,:,5))' ....
    squeeze(xInit_0(2,:,6))' squeeze(xInit_0(2,:,6))' ....
    ];
PolicyRulesStore=vertcat(PolicyRulesStore1,PolicyRulesStore2);
% Iterate on the value function
%load([datapath  'c100.mat'])

for iter=1:Niter
    tic
    disp('Starting Iteration No - ')
    disp(iter)
    IndxSolved=[];
    IndxUnSolved=[];
    
    parfor ctr=1:GridSize
    
    %for ctr=1:GridSize
        u2btild=u2btild_slice(ctr) ;
        R=R_slice(ctr) ;
        s_=s_slice(ctr);
        % At the first pass we use initial guess as the policy functions
        % and the multipliers corresponding to the previous value of
        % coeffecients
        if (ctr==1 || ctr==GridSize)
            xInit=PolicyRulesWorked;
            
        else
            xInit=[PolicyRulesStore(ctr,:)];
        end
        
        [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,xInit',Para);
        
        ExitFlag(ctr)=exitflag;
        VNew(ctr)=V_new;
        PolicyRulesStore(ctr,:)=PolicyRules;
        
        
        
    end
    sprintf(' Done with the parallel computations...it took %1.2f',toc)
    sprintf(' %1.2f  Unresolved so far....',length(find(~(ExitFlag==1))))
    
    
    
    IndxUnSolved=find(~(ExitFlag==1));
    IndxSolved=find(ExitFlag==1);
    IndxSolved_1=IndxSolved(IndxSolved<=GridSize/sSize);
    IndxSolved_2=IndxSolved(IndxSolved>GridSize/sSize);
    
    % Obtain the new coeffecins by projecting the Cheb polynomials for
    % both the value functions
    
    
    cNew(1,:)=funfitxy(V(1),x_state(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
    
    
    cNew(2,:)=funfitxy(V(2),x_state(IndxSolved_2,1:2),VNew(IndxSolved_2)' );
    
    % Store the difference
    cdiff(iter,:)=sum(abs(c-cNew))';
    % update the guess by taking a weighted average of the old and new
    % coeffecients
    c=grelax*cNew+(1-grelax)*c;
    
    disp('Completed Iteration No - ')
    disp(iter)
    
    toc
    sprintf(' %1.0f  Unresolved Points',length(IndxUnSolved))
    
    %PlotFlagPoints
    save([ datapath 'c' num2str(iter)] , 'c', 'cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','x_state','Para','V');
    
end









