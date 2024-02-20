function [x,frac_y,frac_x11]= Within_generation_dynamics_sex_occurs_once(r,f,xtot0,T,muM,k,alpha,Pr,fgenotypes0,tMgenotypes)

% fgenotypes0 is a 1x3 vector embedded in a cell array, defined by {[f00 f01 f10]}.

telapsed=0;  % telapsed allows us to keep track of the current time so we know who should or shouldn't be in meiosis.
tMgenotypesORDERED=sort(tMgenotypes); % tMgenotypesORDERED=sort(tMgenotypes) puts the meiosis time of all strains into ascending order.
S=length(tMgenotypes);  % S is the number of strains in the population.
X=[]; % vector X is used to concatenate the initial population vectors of all strains.

for j=1:S
Xi= xtot0*[fgenotypes0{j} 0 0 0 0 0]*f(j);
X=[X,Xi];
end

for i=1:S+1

  if i==1    
  tspan=[0,tMgenotypesORDERED(i)];    % tspan is the time over which we're running ode45 for each time interval we're concerning e.g.the interval in which m out of S strains are in meiosis
  elseif i==S+1
  tspan=[tMgenotypesORDERED(i-1),T];     
  else
  tspan=[tMgenotypesORDERED(i-1),tMgenotypesORDERED(i)];    
  end

  genotypes_in_meiosis=tMgenotypes(tMgenotypes<=telapsed);  
  genotypes_not_in_meiosis=tMgenotypes(tMgenotypes>telapsed);
  [~,idx_genotypes_in_meiosis] = intersect(tMgenotypes,genotypes_in_meiosis,'stable'); % gives the index of the genotypes that're in meiosis
  [~,idx_genotypes_not_in_meiosis] = intersect(tMgenotypes,genotypes_not_in_meiosis,'stable'); % gives the index of the genotypes that aren't in meiosis



  if S>=2   %
  if tMgenotypes(2)==tMgenotypes(1) && isempty(idx_genotypes_in_meiosis) %
      idx_genotypes_not_in_meiosis=[1,2]; %
  elseif tMgenotypes(2)==tMgenotypes(1) && isempty(idx_genotypes_not_in_meiosis) %
      idx_genotypes_in_meiosis=[1,2]; %
  end %
  end %
  
  
  if tspan(2)==tspan(1)
  else
     [t, x] = ode45(@(t, x) ODE_sex_occurs_once(t,x,r,muM,alpha,k,idx_genotypes_in_meiosis,idx_genotypes_not_in_meiosis,S,Pr), tspan, X);

  end

  X=x(end,:);


  if i==S+1
  telapsed=T;   
  k=1:S;                                                                                        % this 2 lines is the newly added code,
  frac_y=( x(end,k+7*(k-1)+4)+x(end,k+7*(k-1)+5)+x(end,k+7*(k-1)+6)+x(end,k+7*(k-1)+7) );      % the fraction that had sex at time T
  frac_y=sum(frac_y)/sum(x(end,:));
  else
  telapsed=tMgenotypesORDERED(i);
  end


  if i==S+1
  frac_x11=sum((x(end,k+7*(k-1)+3)+x(end,k+7*(k-1)+7)))/sum(x(end,:));
  end
    
end