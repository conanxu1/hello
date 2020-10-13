function []=panding( Ai,taui,n  )

%Ai 以元胞数组形式
%An+1 为A_0

%自适应


	beta=0;
	for i=1:(n+1)
		beta=beta+norm(Ai{i,1});
	end

	[dimxm,dimxn]=size(Ai{i,1});

	h0=1e-5;
	pz=zeros(1,round(1e8));
	
	i=1;
	taunow=	 0;
	
	while(    taunow<0.5    )
		
		snow=	beta*cos( (taunow-1/4)*2*pi ) +1i*beta*sin( (taunow-1/4)*2*pi) ;

					tp=snow*eye(dimxm);
					tp=tp-Ai{n+1,1};
					for j=1:n
						tp=tp- Ai{j,1}*exp(-snow *taui(j)  );
					end
					pz(i)=det(tp);
				
				
										if norm(pz(i) )>0.1
											h=h0;
										else 
											if norm(pz(i)) >0.01
												h=h0*1e-2;
											else  
												if norm(pz(i)) >1e-9
													h=h0*1e-4;
												else
													if norm(pz(i)) >1e-11
														h=h0*1e-6;
													else
														norm(pz(i))
														fprintf('no')
														return 
													end
												end
											end
										end
	taunow=taunow+h;
	i=i+1;
	end
	
	
	while(    taunow<=1    )
		
		snow=	 0+1i*  (beta-  2*beta*   (taunow-1/2)*2 );

					tp=snow*eye(dimxm);
					tp=tp-Ai{n+1,1};
					for j=1:n
						tp=tp- Ai{j,1}*exp(-snow *taui(j)  );
					end
					pz(i)=det(tp);
				
										if norm(pz(i) )>0.1
											h=h0;
										else 
											if norm(pz(i)) >0.01
												h=h0*1e-2;
											else  
												if norm(pz(i)) >1e-9
													h=h0*1e-4;
												else
													if norm(pz(i)) >1e-11
														h=h0*1e-6;
													else
													
														norm(pz(i))
														fprintf('no')
														return 
													end
												end
											end
										end
	taunow=taunow+h;
	i=i+1;
	end
	
	pzend=pz(1,1:(i-1));
	
	

	
	Argpz=angle(pzend);
	Num=0;
	for t=2:length(Argpz)
		if abs(Argpz(t)-Argpz(t-1))>1
			if (Argpz(t)-Argpz(t-1))>0
				Num=Num-1;
				Argpz(t:end)=Argpz(t:end)-2*pi;
			else
				Num=Num+1;
				Argpz(t:end)=Argpz(t:end)+2*pi;
			end
		end
	end
	
	W=  round(   abs(Argpz(length(Argpz))-Argpz(1))/(2*pi)  )
	
	
	% plot(Argpz)
	% hold on
	% plot(angle(pzend),'*');
	figure(109)
	plot(pzend,'*')
	                
	
	
	
	
	
	
	
	
	
	
	
	

	% % % % % % % length(pzend)
	% % % % % % % normpz=pzend;
    % % % % % % % for i=1:length(normpz)      
      % % % % % % % normpz(i)=norm(pzend(i));
    % % % % % % % end

    % % % % % % % min(normpz)




    
    
    
    
	% % % % % h=0.00001;
	% % % % % tau1=0:h:1/2;
	% % % % % tau2=(1/2+h):h:1;
	% % % % % y1=  beta*cos( (tau1-1/4)*2*pi ) +1i*beta*sin( (tau1-1/4)*2*pi)                ;
	% % % % % y2=    0+1i*  (beta-( (tau2-1/2)*2 ) *2*beta);
	% % % % % ss=[y1 ,y2];
	
	
	
	
	
	
	% % % % % outp=zeros(1,length(ss));
	% % % % % [dimxm,dimxn]=size(Ai{1,1});
	% % % % % for i=1:length(ss)
		% % % % % tp=ss(i)*zeros(dimxm);
		% % % % % tp=tp-Ai{n+1,1};
		% % % % % for j=1:n
			% % % % % tp=tp- Ai{j,1}*exp(-ss(i) *taui(j)  );
		% % % % % end
		% % % % % outp(i)=det(tp);
    % % % % % end
	
    
    
    
    
% % % % %     plot(outp,'*')



     % % % % % normop=outp;




    % % % % % for i=1:length(outp)      
      % % % % % normop(i)=norm(outp(i));
    % % % % % end

    % % % % % min(normop)

    % % % % % plot(outp,'*')

% % % % % for i=1:length(outp)
% % % % % 	plot(outp(i),'*')
% % % % % %     pause(0.001);
% % % % %     hold on
% % % % %     end







end
