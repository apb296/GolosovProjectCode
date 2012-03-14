
%BuildGrid

% 1. u2*tilde_b_2
u2btildMin=1;
u2btildMax=200;

% Use this code to get a grid that is dense near the end points
% BrPoint1=u2btildMin+(u2btildMax-u2btildMin)/5;
% BrPoint2=u2btildMin+(u2btildMax-u2btildMin)*4/5;
% if (u2btildGridSize/5)<2
% u2btildGrid=[linspace(BrPoint1*.98,BrPoint2,4*u2btildGridSize/5)];
% u2btildGrid= [u2btildMin u2btildGrid u2btildMax];
% else
% u2btildGrid=linspace(u2btildMin,BrPoint1,round(u2btildGridSize/5));
% u2btildGrid=[u2btildGrid linspace(BrPoint1*.98,BrPoint2,4*u2btildGridSize/5)];
% u2btildGrid= [u2btildGrid linspace(BrPoint2,u2btildMax,round(u2btildGridSize/5))];
% end
% 
% u2btildGridSize=length(u2btildGrid);
% Para.u2btildGrid=u2btildGrid;
% Para.u2btildGridSize=u2btildGridSize;
u2btildGrid=linspace(u2btildMin,u2btildMax,u2btildGridSize);
Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=-2;
Para.u2btildUL=u2btildMax*1.05;

% R=u_2/u_1 = (c1/c2)^(sigma)
RMin=4;
RMax=8;
RGrid=linspace(RMin,RMax,RGridSize);
Para.RGrif=RGrid;
GridSize=u2btildGridSize*RGridSize*sSize;
Para.GridSize=GridSize;

%% Define the funtional space

V(1) = fundefn('cheb',[OrderOfAppx_u2btild OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
V(2) = V(1);

GridPoints=[u2btildMin u2btildMax;RMin RMax];
rowLabels = {'$x$','$R$'};
columnLabels = {'Lower Bound','Upper Bounds'};
matrix2latex(GridPoints, [texpath 'GridPoints.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
