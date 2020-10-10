close all
clear all
clc
global A
global B
global Atau

R1=553.1;
R2=1069.8;


C1=0.670*1e-3;
C2=0.620*1e-3;
Rd1=553.1;
Rd2=553.1;


A=   -[1/C1/R1  0;0   1/C2/R2];


Atau=  -[1/C1/Rd1    ,  0;
				1/C2/Rd1    , 1/C2/Rd2];

B=[1/C1;  1/C2]/1000;
 
tau=0.2 ;
tspan=10;
 
Delta= 0.01;
numtau=round(tau/Delta);

global ttt

global canshu
tt=0:0.01:tspan;
ttt=tt;
canshu=rand(1,length(tt));


 
options = ddeset('InitialStep',1e-13,'Jumps',0,'AbsTol',1e-13,'InitialY',[0;0],'MaxStep' , 2*1e-3  ); 
sol = dde23(@rea,[tau],@reahis,[0, tspan],options);    %%delay  history span


 
tt=tt.';
xx=deval(sol,tt).';


uu=zeros(length(tt),1);
 for i=1:length(tt)
 uu(i,1)=myut(tt(i));
 end



plot(tt,uu,'o-');
hold on
plot(tt,xx(:,1)+xx(:,2),'+-');
hold on
plot(tt,xx(:,2),'*-');
hold on






% % % y1=zeros(   (length(xx(:,1)) -numtau)  ,1);
% % % for i=(numtau+1):length(xx(:,1))
	% % % y1(i-numtau,1)=xx(i,1)+xx(i,2)-xx(i-numtau,1);
% % % end

% % % y2=zeros(   (length(xx(:,1)) -numtau)  ,1);
% % % for i=(numtau+1):length(xx(:,1))
	% % % y2(i-numtau,1)=xx(i,2)-xx(i-numtau,2);
% % % end



% % % plot(  ((numtau+1):length(xx(:,1)))*Delta,y1,'og-'             );
% % % hold on
% % % plot(  ((numtau+1):length(xx(:,1)))*Delta,y2,'xr-'             );
% % % hold on


% % % plot(  [0,tspan],[0,0],'x-'             );








 
 % [tt,yy] = ode45(@ff, [0,5], [0;0]);
 
 
 
 
 
 % uu=zeros( length(tt) , 1  );
 
 % for i=1:length(tt)
 % uu(i,1)=myut(tt(i));
 % end
 
 % plot(tt, uu, '+-');
 % hold on

% plot(tt, yy(:,1)+yy(:,2), 'o-');
% hold on

% plot(tt, yy(:,1), '*-');
 