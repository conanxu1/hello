function dydt =rea(t,y,Z)
global A
global Atau
global B 

 
 
dydt =A*y+Atau*Z(:,1)   +B*myut(t);
 
end