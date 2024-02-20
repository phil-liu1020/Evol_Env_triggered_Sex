function [tM,frac_had_sex] = Evolutionary_Dynamics_sex_occurs_once(tM0, T, r,alpha,f0,xtot0,mu,muM,deltatM, NEVOL,k,Pr)

% fgenotypes0 is the initial frequency of each of the x00 and x1 genotypes in the simulation. The input should be e.g fgenotypes0={[f00 f1]}. Only a 1D array please.

% tM0 - initial trait value of the meiosis time

% f0 - initial frequency of mutant with a different tM

% muM - rate of genotype mutation into one another i.e x00 into x1

% mu - rate at which mutant traits with different tMs emerge

% xtot0 - initial total population min the system. This includes all trait
% values.

% Ngens is assumed to be 300 here, the number of generations to run the
% equilibration for.

%frac_had_sex  is the fraction of individuals that had the opportunity to
%engage in sex. This quantity is only checked at the NEVOLth generation.

tM=zeros(1,NEVOL);
fgenotypes=equilibriate_frequency_of_genotypes_sex_occurs_once(500,r,T,xtot0,mu,k,tM0,alpha,Pr);
fgenotypes0={[fgenotypes(1,end),fgenotypes(2,end)/2,fgenotypes(2,end)/2]};

tMgenotypesDATA=cell(2,NEVOL);
tMgenotypes=tM0;
Nev=1;
f=1;
fgenotypesARRAYend={fgenotypes0{1,1}};
f00DATA=cell(1,NEVOL);
fraction_lethal_genotypes=zeros(1,NEVOL);

while Nev<=NEVOL

   
       keepstrains=0;    
       rr=rand;
       prob = f;
       if any(f<f0)
             index_to_exclude = f<f0; % 
             logical_index = true(size(f)); % 
             logical_index(index_to_exclude) = false; % 
             new_vector = f(logical_index); % 
             index_list=(1:1:length(f))';
             index_list=index_list.*logical_index;
             index_list(index_list==0)=[];
             q= sum(rr >= cumsum([0;new_vector/sum(new_vector)])); %fine here, no error
             q = index_list(q);
       else
       q= sum(rr >= cumsum([0;prob/sum(f)]));    % q is the index of 
       end


       tMgenotypes(length(tMgenotypes)+1)=round(tMgenotypes(q)+deltatM*sign(randn),2);


       if any(tMgenotypes>T)
       tMgenotypes(tMgenotypes>T)=T;
       end


       S=length(tMgenotypes);
       fgenotypestoconcatenate=fgenotypesARRAYend{1,q};

       if ~any(groupcounts(tMgenotypes')>1)                  % Here also REVISED
       fgenotypesARRAYend{1,S}=fgenotypestoconcatenate;
       end

           % make sure tMgenotypes is a column vector
           if iscolumn(tMgenotypes)
               tMgenotypes=tMgenotypes';
           end
           
           if any(groupcounts(tMgenotypes')>1)                                   % This if statement is to merge mutant strains with same mass as an existing strain.
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
       
       G=geornd(mu)+1;

       if Nev<NEVOL
       [f,fgenotypesARRAYend,~]=Invasion_dynamics_sex_occurs_once(tMgenotypes, T, r, alpha, f,xtot0,G,k,muM,Pr,fgenotypesARRAYend);
       elseif Nev==NEVOL
       [f,fgenotypesARRAYend,frac_had_sex]=Invasion_dynamics_sex_occurs_once(tMgenotypes, T, r, alpha, f,xtot0,G,k,muM,Pr,fgenotypesARRAYend);    
       end

       f00DATA{Nev}=cellfun(@(x) x(1), fgenotypesARRAYend(1:S));



       keepstrains=zeros(S,1);
       keepstrains=(f>=10^(-2));

f(keepstrains==0)=[];
f=f/sum(f);
tMgenotypes([keepstrains]==0)=[];
tMgenotypesDATA{1,Nev}=tMgenotypes;
tMgenotypesDATA{2,Nev}=f;

if length(f)==1
    f=1;
end 

       tM(Nev)=sum(tMgenotypes'.*f);


fprintf('Processing %d...',Nev);

Nev=Nev+1;
end