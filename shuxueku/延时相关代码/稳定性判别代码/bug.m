close all
clear all
clc

N=1
dimx=2
 
 
A_i=cell(N+1,1);
A_i{1,1}= -ones(dimx,dimx)+triu(sin(ones(dimx,dimx)))- tril(sin(ones(dimx,dimx)))    ;
A_i{2,1}= -0.2*triu(sin(ones(dimx,dimx)))  -9*ones(dimx,dimx)  ;





tau_i(1)=0.8;



panding(A_i,tau_i,N)



