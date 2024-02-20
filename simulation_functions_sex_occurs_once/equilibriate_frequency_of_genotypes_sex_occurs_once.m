function f= equilibriate_frequency_of_genotypes_sex_occurs_once(Ngens,r,T,xtot0,mu,k,tM,alpha,Pr)

if tM==T
x = Within_generation_dynamics_sex_occurs_once(r,1,xtot0,T,mu,k,alpha,Pr,{[1 0 0]},T*0.9999);
else
x = Within_generation_dynamics_sex_occurs_once(r,1,xtot0,T,mu,k,alpha,Pr,{[1 0 0]},tM);
end


f(1,1)=( x(end,1)+x(end,5) )/( x(end,1)+x(end,2)+x(end,3)+x(end,5)+x(end,6)+x(end,7) );
f(2,1)=( x(end,2)+x(end,3)+x(end,6)+x(end,7) )/( x(end,1)+x(end,2)+x(end,3)+x(end,5)+x(end,6)+x(end,7) );

for i=2:Ngens

if tM==T
x = Within_generation_dynamics_sex_occurs_once(r,1,xtot0,T,mu,k,alpha,Pr,{[f(1,i-1) f(2,i-1)/2 f(2,i-1)/2]},T*0.9999);
else  
x = Within_generation_dynamics_sex_occurs_once(r,1,xtot0,T,mu,k,alpha,Pr,{[f(1,i-1) f(2,i-1)/2 f(2,i-1)/2]},tM);
end 

f(1,i)=( x(end,1)+x(end,5) )/( x(end,1)+x(end,2)+x(end,3)+x(end,5)+x(end,6)+x(end,7) );
f(2,i)=( x(end,2)+x(end,3)+x(end,6)+x(end,7) )/( x(end,1)+x(end,2)+x(end,3)+x(end,5)+x(end,6)+x(end,7) );
end