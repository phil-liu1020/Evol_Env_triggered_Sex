function [gend,fgenotypesARRAYend] = Invasion_dynamics_sex_arbitrary(tMgenotypes, T, r0,r1,r2, alpha, f,xtot0,G,k,muM,Pc,fgenotypes0)

S=length(tMgenotypes);
g=zeros(S,1);
fgenotypesARRAY=cell(G,1);
fgenotypesARRAY{1}=fgenotypes0;
 
for i=1:G
    
   x = Within_generation_Dynamics_sex_arbitrary(r0,r1,r2, f,xtot0,muM,k,alpha,tMgenotypes,T,Pc,fgenotypesARRAY{i});
   
   [f,fgenotypesARRAY{i+1}]= Survival_function_sex_arbitrary(x);
   
   g(:,i)=f;

end

gend=g(:,end);  
fgenotypesARRAYend=fgenotypesARRAY{end};