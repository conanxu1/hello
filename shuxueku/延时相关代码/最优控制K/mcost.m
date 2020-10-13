function val=mcost(canshu)
[m,n]=size(canshu);
global dimx
global dimu
global jieshu
global x0k
global wwk
global gsool
global Ai
global B


global tauN
global tspan

if n>1
canshu=canshu.';
end


KK=reshape(canshu,dimu,dimx*(tauN+1) );
for i=1:(tauN+1)
	Ki{i,1}=KK(:,((i-1)*dimx+1) : ((i-1)*dimx+dimx)     );
end

jijiezhen(Ki)


Ftk=zeros(dimx,dimx);


val=0;
for i=1:(jieshu+1)

	for j=1:dimx
		Ftk(:,j)=deval(gsool{j,1} ,x0k(i));
	end

	val=val+wwk(i)*norm( Ftk)*expm( x0k(i));
end

end