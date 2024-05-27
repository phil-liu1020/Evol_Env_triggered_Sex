function [tM,tMgenotypesDATA] = Evolutionary_Dynamics_sex_arbitrary(tM0, T, r0,r1,r2,alpha,f0,xtot0,mu,muM,deltatM, NEVOL,k,Pr)

% fgenotypes0 is the initial frequency of each of the x00, x01 and x10
% genotypes in the simulation. The input should be e.g fgenotypes0={[f00
% f01 f10]}. Only a 1D array please.

% tM0 - initial trait value of the meiosis time

% f0 - initial frequency of mutant with a different tM

% muM - rate of genotype mutation into one another i.e x00 into x01 and x10

% mu - rate at which mutant traits with different tMs emerge

% xtot0 - initial total population min the system. This includes all trait
% values.

tM=zeros(1,NEVOL);
fgenotypes=equilibriate_frequency_of_genotypes_sex_arbitrary(1000,r0,r1,r2,T,xtot0,mu,k,tM0,alpha,Pr);
fgenotypes0={[fgenotypes(1,end),fgenotypes(2,end)/2,fgenotypes(2,end)/2]};

tMgenotypesDATA=cell(2,NEVOL);
tMgenotypes=tM0;
Nev=1;
f=1;
tMcur=0;
fgenotypesARRAYend={fgenotypes0{1,1}};

while Nev<=NEVOL
   
       keepstrains=0;    
       rr=rand;
       prob = f;
       if any(f<f0)
             index_to_exclude = f<f0;  
             logical_index = true(size(f));  
             logical_index(index_to_exclude) = false;
             new_vector = f(logical_index); 
             index_list=(1:1:length(f))';
             index_list=index_list.*logical_index;
             index_list(index_list==0)=[];
             q= sum(rr >= cumsum([0;new_vector/sum(new_vector)])); 
             q = index_list(q);
       else
       q= sum(rr >= cumsum([0;prob/sum(f)]));                              
       end


       tMgenotypes(length(tMgenotypes)+1)=round(tMgenotypes(q)+deltatM*sign(randn),1);

       if any(tMgenotypes>T)
       tMgenotypes(tMgenotypes>T)=T;
       end


       S=length(tMgenotypes);
       fgenotypestoconcatenate=fgenotypesARRAYend{1,q};

       if ~any(groupcounts(tMgenotypes')>1)                
       fgenotypesARRAYend{1,S}=fgenotypestoconcatenate;
       end

           % make sure tMgenotypes is a column vector
           if iscolumn(tMgenotypes)
               tMgenotypes=tMgenotypes';
           end
           
           if any(groupcounts(tMgenotypes')>1)                                   % This if statement is to merge mutant strains with same tM as an existing strain.
           u=unique(tMgenotypes);
           u=u(1<histc(tMgenotypes,unique(tMgenotypes))); 
           fgenotypeofDUPLICATE=fgenotypesARRAYend{:,find(tMgenotypes==u, 1 )}*f(find(tMgenotypes==u, 1 ),1)+fgenotypestoconcatenate*f0;  %REVISED piece of code
           fgenotypeofDUPLICATE=fgenotypeofDUPLICATE/sum(fgenotypeofDUPLICATE);
           
           f(find(tMgenotypes==u, 1 ),1)=f(find(tMgenotypes==u, 1 ),1)+f0;         
           f(q)=f(q)-f0;                                                         % largest index of the element that occurs more than once
           tMgenotypes(find(tMgenotypes==u, 1, 'last' ) )=[];
           fgenotypesARRAYend(find(tMgenotypes==u, 1 )) = {fgenotypeofDUPLICATE};      %This piece of code is the REVISED piece to concatenate frequencies.
           S=length(tMgenotypes);
           else
           f(end+1,1) = f0;
           f(q)=f(q)-f0;
           end 

       %-----------------------------------------------------
       if any(f<0) %%%% Here's a code for a sanity check!!!
    keyboard       %%%% Here's a code for a sanity check!!!
       end         %%%% Here's a code for a sanity check!!!
       %-----------------------------------------------------

       G=geornd(mu)+1;
       [f,fgenotypesARRAYend]=Invasion_dynamics_sex_arbitrary(tMgenotypes, T, r0,r1,r2, alpha, f,xtot0,G,k,muM,Pr,fgenotypesARRAYend);

       %-----------------------------------------------------
       if any(f<0) %%%% Here's a code for a sanity check!!!
    keyboard       %%%% Here's a code for a sanity check!!!
       end         %%%% Here's a code for a sanity check!!!
       %-----------------------------------------------------

       keepstrains=zeros(S,1);
       for i=1:length(f)
           if f(i)<10^(-2)
           keepstrains(i,:)=0;
           else
           keepstrains(i,:)=1;    
           end
       end



f(keepstrains==0)=[];
f=f/sum(f);
tMgenotypes([keepstrains]==0)=[];
tMgenotypesDATA{1,Nev}=tMgenotypes;
tMgenotypesDATA{2,Nev}=f;

if length(f)==1
    f=1;
end
       

       for i=1:length(f)        % mean mass of the population
       tMcur=tMcur+f(i)*tMgenotypes(i);
       end       
       tM(Nev)=tMcur;
       tMcur=0;

fprintf('Processing %d...',Nev);

Nev=Nev+1;
end