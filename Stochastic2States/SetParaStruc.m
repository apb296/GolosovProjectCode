% Set Params


if strcmp(computer,'PCWIN')
    sl='\';
    coresize=2;
else
    sl='/';
    coresize=8;
end


    
% 1. Paramters describing the preferences
theta_1=1; % type of Agent 1
theta_2=0; % Type of Agent 2
sigma=1; % Risk aversion
sigma_1=sigma;
sigma_2=sigma;
gamma=2; % Inv Frish Elasticity of labor
beta=.98 ;% subjective time discount factor;

% 2. Technology
g_l=.18; % Government expenditure in low state s_l
g_h=.2; % Government expenditure in high state s_h
P=[.5 .5;.5 .5]; % Transition Matrix for g shocks
alpha_1=.5;
alpha_2=1-alpha_1;
alpha=[alpha_1 alpha_2]; % Pareto Weights for agent 1 and Agent 2;
sSize=2; % Dimension of the markov state
% 3. Others
pertub=0.00;
ctol=1e-7;
grelax=.9;
Niter=500;

  u2btildGridSize=20;
  RGridSize=8;
  
  OrderOfAppx_u2btild=4;
  OrderOfApprx_R=3;
 
compeconpath=[pwd sl 'compecon2011' sl];
knitropath=[pwd sl 'knitro' sl];
texpath= [pwd sl 'Tex' sl] ;
plotpath= [pwd sl 'Graphs' sl] ;
datapath=[pwd sl 'Data' sl] ;
mkdir(texpath)
mkdir(plotpath)
mkdir(datapath)
addpath(genpath(compeconpath))
addpath(genpath(knitropath))


rowLabels = {'$\sigma$','$\gamma$','$\beta$', '$g_{l}$','$g_{h}$'};
columnLabels = {};
matrix2latex([sigma gamma beta g_l g_h]', [texpath 'Param.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
rowLabels = {'$\theta$','$\alpha$'};
columnLabels = {'Agent 1','Agent 2'};
matrix2latex([theta_1 alpha_1 ;theta_2 alpha_2]', [texpath 'AgentParam.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
rowLabels={'$g_l$','$g_h$'};
columnLabels={'$g_l$','$g_h$'};
matrix2latex(P, [texpath 'Pr.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
Para.ctol=ctol;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.sigma=sigma;
Para.sigma_1=sigma_1;
Para.sigma_2=sigma_2;
Para.gamma=gamma;
Para.beta=beta ;
Para.g_l=g_l;
Para.g_h=g_h;
g=[g_l g_h];
Para.g=g;
Para.P=P;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
