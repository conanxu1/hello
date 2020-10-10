close all
clear all
clc

beta=2;

	h=0.001;
	
	tau1=0:h:1/2;
	tau2=(1/2+h):h:1;
	y1=  beta*cos( (tau1-1/4)*2*pi ) +1i*beta*sin( (tau1-1/4)*2*pi)                ;
	y2=    0+1i*  (beta-( (tau2-1/2)*2 ) *2*beta);
	
	ss=[y1 ,y2];
	
	for i=1:length(y1)
		r(i)=norm(y1(i));
	end	
	plot(y1);
	hold on
	plot(y2);
	
	max(diff([tau1,tau2]))