function [] =jijiezhen(Ki)
global dimx
global tspan
global Ft
global tauk
global temAk
global gsool
global tauN
global Ai
global B



global ee1
global ee2
global ee3


for i=1:(tauN+1)
	temAk{i,1}=Ai{i,1}+B*Ki{i,1};
end

gsool=cell(dimx,1);


% global ttdelta


for i=1:dimx

	x0=zeros(dimx,1);
	x0(i)=1;
	options = ddeset('InitialStep',ee1,'Jumps',0,'AbsTol',ee2,'InitialY',x0,'MaxStep',ee3 ); 
	
	
	gsool{i,1} = dde23(@solFt,tauk,@solFthis,[0, tspan],options);    %%delay  history span

% yy = deval(sol,ttdelta);
% Ft(:,i,:)=yy;



end





end

