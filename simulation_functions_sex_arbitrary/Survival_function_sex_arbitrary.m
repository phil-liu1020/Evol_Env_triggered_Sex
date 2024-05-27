function [f,fgenotypes]= Survival_function_sex_arbitrary(x)

S=length(x(end,:))/4;

fgenotypes=cell(1,S);

for i=1:S
   a(i)= x(end,i+3*(i-1)) + x(end,i+3*(i-1)+1) + x(end,i+3*(i-1)+2);
   fgenotypes{i}=[x(end,i+3*(i-1))/a(i) x(end,i+3*(i-1)+1)/a(i) x(end,i+3*(i-1)+2)/a(i)];
end


for i=1:S
    f(i)=(   x(end,i+3*(i-1)) + x(end,i+3*(i-1)+1) + x(end,i+3*(i-1)+2)  )/sum(a);
end