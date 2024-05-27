function x = Within_generation_Dynamics_sex_arbitrary(r0,r1,r2, f,xtot0,muM,k,alpha,tMgenotypes,T,Pr,fgenotypes0) %, tMgenotypes X and Y are vectors 

telapsed=0;
tMgenotypesORDERED=sort(tMgenotypes);
S=length(tMgenotypes);
X=[];


for j=1:S
Xi= xtot0*[fgenotypes0{j} 0]*f(j);
X=[X,Xi];
end

  for i=1:S+1
  
  if i==1    
  tspan=[0,tMgenotypesORDERED(i)];    % time over which we're running ode45
  elseif i==S+1
  tspan=[tMgenotypesORDERED(i-1),T];     
  else
  tspan=[tMgenotypesORDERED(i-1),tMgenotypesORDERED(i)];    
  end
  
  genotypes_in_meiosis=tMgenotypes(tMgenotypes<=telapsed);
  genotypes_not_in_meiosis=tMgenotypes(tMgenotypes>telapsed);
  [~,idx_genotypes_in_meiosis] = intersect(tMgenotypes,genotypes_in_meiosis,'stable');
  [~,idx_genotypes_not_in_meiosis] = intersect(tMgenotypes,genotypes_not_in_meiosis,'stable');
  


  if S>=2   
  if tMgenotypes(2)==tMgenotypes(1) && isempty(idx_genotypes_in_meiosis) 
      idx_genotypes_not_in_meiosis=[1,2]; 
  elseif tMgenotypes(2)==tMgenotypes(1) && isempty(idx_genotypes_not_in_meiosis) 
      idx_genotypes_in_meiosis=[1,2]; 
  end 
  end 
  
  
  if tspan(2)==tspan(1)
  else
      [t, x] = ode45(@(t, x)ODE_sex_arbitrary(t,x,r0,r1,r2,muM,alpha,k,idx_genotypes_in_meiosis,idx_genotypes_not_in_meiosis,S,Pr), tspan, X);
  end





  
  X=x(end,:);
  
  if i==S+1
  telapsed=T;    
  else
  telapsed=tMgenotypesORDERED(i);
  end
    
  end