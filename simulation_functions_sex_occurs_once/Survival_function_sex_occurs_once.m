function [f,fgenotypes]= Survival_function_sex_occurs_once(x)

% x(end,i+7*(i-1)) is x00      % x(end,i+7*(i-1)+3) is x11     % x(end,i+7*(i-1)+6) is y10
% x(end,i+7*(i-1)+1) is x01     % x(end,i+7*(i-1)+4) is y00    % x(end,i+7*(i-1)+7) is y11
% x(end,i+7*(i-1)+2) is x10     % x(end,i+7*(i-1)+5) is y01

S=length(x(end,:))/8; 

fgenotypes=cell(1,S);

for i=1:S
   a(i)= x(end,i+7*(i-1)) + x(end,i+7*(i-1)+1) + x(end,i+7*(i-1)+2) + x(end,i+7*(i-1)+4) + x(end,i+7*(i-1)+5) + x(end,i+7*(i-1)+6);
   fgenotypes{i}=[(x(end,i+7*(i-1))+x(end,i+7*(i-1)+4))     (x(end,i+7*(i-1)+1)+x(end,i+7*(i-1)+5))     (x(end,i+7*(i-1)+2)+x(end,i+7*(i-1)+6))]/a(i);
end


for i=1:S
    f(i)=(  x(end,i+7*(i-1)) + x(end,i+7*(i-1)+1) + x(end,i+7*(i-1)+2) + x(end,i+7*(i-1)+4) + x(end,i+7*(i-1)+5) + x(end,i+7*(i-1)+6)  )/sum(a);
end        