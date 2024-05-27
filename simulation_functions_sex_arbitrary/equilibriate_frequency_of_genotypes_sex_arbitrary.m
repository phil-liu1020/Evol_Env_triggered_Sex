function f= equilibriate_frequency_of_genotypes_sex_arbitrary(Ngens,r0,r1,r2,T,xtot0,mu,k,tM,alpha,Pr)

if tM==T
x = Within_generation_Dynamics_sex_arbitrary(r0,r1,r2, 1,xtot0,mu,k,alpha,T*0.9999,T,Pr,{[1 0 0]});
else
x = Within_generation_Dynamics_sex_arbitrary(r0,r1,r2, 1,xtot0,mu,k,alpha,tM,T,Pr,{[1 0 0]});
end


f(1,1)=x(end,1)/sum(x(end,(1:3)));
f(2,1)=(x(end,2)+x(end,3))/sum(x(end,(1:3)));

for i=2:Ngens

if tM==T
x = Within_generation_Dynamics_sex_arbitrary(r0,r1,r2, 1,xtot0,mu,k,alpha,T*0.9999,T,Pr,{[f(1,i-1) f(2,i-1)/2 f(2,i-1)/2]});
else
x = Within_generation_Dynamics_sex_arbitrary(r0,r1,r2, 1,xtot0,mu,k,alpha,tM,T,Pr,{[f(1,i-1) f(2,i-1)/2 f(2,i-1)/2]});
end

f(1,i)=x(end,1)/sum(x(end,(1:3)));
f(2,i)=(x(end,2)+x(end,3))/sum(x(end,(1:3)));
end