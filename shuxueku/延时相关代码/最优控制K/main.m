close all
clear all
clc
global x0k
global wwk


global jieshu
global dimx
global dimu
global tauN

global Ai
global B
global tspan
global tauk
 

global ee1
global ee2
global ee3

ee1=1e-8;
ee2=1e-8;
ee3=8*1e-2;


dimx=2;
dimu=1;
tauN=1;

tauk=[1];

jieshu=11;
[x0k,wwk]=xishu(jieshu);



tspan=x0k(end,1)


%%%%%%%%%%%%%%%


R1=553.1;
R2=1069.8;


C1=0.670*1e-3;
C2=0.620*1e-3;
Rd1=553.1;
Rd2=553.1;

Ai=cell(tauN+1,1);
A=   -[1/C1/R1  0;0   1/C2/R2];


Atau=  -[1/C1/Rd1    ,  0;
				1/C2/Rd1    , 1/C2/Rd2];

B=[1/C1;  1/C2]/1000;


Ai{1,1}=Atau;
Ai{2,1}=A;


nnum=dimu*dimx*(tauN+1);
LB=-10*ones(nnum,1);
UB= 10*ones(nnum,1);


options=gaoptimset('paretoFraction',0.4,'populationsize',60,'generations',3,'stallGenLimit',50,'TolFun',1e-10,'PlotFcns',@gaplotbestf);

[thth,coco]=ga(  @(canshu) mcost(canshu) , nnum  ,[],[],[],[],LB,UB,[],options);

KKi=canshu2Ki(thth);
newAi=cell(tauN+1,1);
for i=1:(tauN+1)
	newAi{i,1}=Ai{i,1}+B*KKi{i,1};	
end
panding(newAi,tauk,tauN);


% Ki=cell(tauN+1,1);


% Ki{1,1}=eye(dimu,dimx);
% Ki{2,1}=-eye(dimu,dimx);

% tic

% toc
