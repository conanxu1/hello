function [xk,Ak]=xishu(n)


pn=peidian(n);
pn2=peidian(n+1);

xk=roots(pn);

xk=sort(xk);


Ak=zeros(1,n+1);

for i=0:n
	Ak(1,i+1)=  factorial(n+1)^2/ ( polyval(pn2,xk(i+1)))^2   * xk(i+1) ;
end

end