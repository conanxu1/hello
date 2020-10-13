function KKi=canshu2Ki(canshu)
global dimx
global dimu
global tauN

KKi=cell(tauN+1,1);
KK=reshape(canshu,dimu,dimx*(tauN+1) );
for i=1:(tauN+1)
	KKi{i,1}=KK(:,((i-1)*dimx+1) : ((i-1)*dimx+dimx)     );
end



end