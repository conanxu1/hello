function dydt = solFt(t,y,Z)
global dimx
global temAk
 
global tauN


    dydt =temAk{tauN+1,1}*y;
    for i=1:tauN
        dydt =dydt+temAk{i,1}*Z(:,i);
    end

end