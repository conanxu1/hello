function y =myut(t)

global canshu
global ttt


zhouqi=1;
tt=t-zhouqi*floor(t/zhouqi) ;

A=0.5;




y=0;
if t<0.2
    
    y=A;
else
    y=0.01*interp1(ttt,canshu,t,'linear')       ;
end
    
% 
% if tt<0.2
% 	y=A;
% 
%     else if tt<=zhouqi
% 	y=-A;
% 	end
% end

 


% y=0.5*sin(2*pi*1*t);


end