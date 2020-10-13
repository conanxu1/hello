 
function pnnext=peidian(n) 

n=n+1;
gpn=zeros(1,n+1);
gpn1=zeros(1,n+1);

pn=[0,-1,1];
gpn(1,(end-2):end)=pn;

pn1=[1,-4,2];
gpn1(1,(end-2):end)=pn1;


for i=2:(n-1)
	% pnnext=
	
	tp=conv(gpn1,[-1,0]);
	tp=tp(1,2:end);
	
	
	pnnext=(2*i+1)*gpn1-i^2*gpn+tp;
	
	
	gpn=gpn1;
	gpn1=pnnext;


end




end

% y=polyder(xishu)