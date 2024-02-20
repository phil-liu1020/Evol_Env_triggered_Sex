function [gend,fgenotypesARRAYend,frac_had_sex] = Invasion_dynamics_sex_occurs_once(tMgenotypes, T, r, alpha, f,xtot0,G,k,muM,Pr,fgenotypes0)

S=length(tMgenotypes);
g=zeros(S,1);    % g is the vector of frequencies of each strain with a different trait value. It's initialised as "zeros(S,1)", but at the end it will be an SxNEVOL vector, where each row shows the frequency of each trait value. 
fgenotypesARRAY=cell(G,1); % the function fgenotypesARRAY records the frequency of f00 and f1 at the end of every generation i.e at every value of t_g (t_g is the generational timescale). Each array element is a 1x2 vector.
fgenotypesARRAY{1}=fgenotypes0; % fgenotypes0 is the initial frequency of the genotypes f00 and f1.
 
for i=1:G
    
   [x,frac_y] = Within_generation_dynamics_sex_occurs_once(r,f,xtot0,T,muM,k,alpha,Pr,fgenotypesARRAY{i},tMgenotypes);
   
   [f,fgenotypesARRAY{i+1}]= Survival_function_sex_occurs_once(x);
   
   g(:,i)=f;
   
end

gend=g(:,end);  %gend is the frequency of each strain right before the next mutant is introduced
fgenotypesARRAYend=fgenotypesARRAY{end}; %gend is the frequency of genotypes f00 and f1 right before the next mutant is introduced
frac_had_sex=frac_y;