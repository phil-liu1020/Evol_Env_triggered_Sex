function dx = ODE_sex_occurs_once(t,x,r,muH,alpha,k,idx_genotypes_in_meiosis,idx_genotypes_not_in_meiosis,S,Pr)

% Here, swapping of the top locus occurs with probability (1-k)Pr,
% swapping of the middle and bottom loci occur each with probability 0.5*k*Pr,


m=2/3;
dx=zeros(8*S,1);

if ~isempty(idx_genotypes_in_meiosis)
    
    for i=1:length(idx_genotypes_in_meiosis)
    d=0; e=0; h=0; 

    a=x(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1));
    b=x(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+1);
    c=x(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+2);
    

    for j=1:length(idx_genotypes_in_meiosis)

         if j~=i
         xx=x(idx_genotypes_in_meiosis(j)+7*(idx_genotypes_in_meiosis(j)-1)+1);
         d=d+xx;
         yy=x(idx_genotypes_in_meiosis(j)+7*(idx_genotypes_in_meiosis(j)-1)+2);
         e=e+yy;    
         zz=x(idx_genotypes_in_meiosis(j)+7*(idx_genotypes_in_meiosis(j)-1));
         h=h+zz;
         end

    end
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1))=-alpha*a^2-alpha*a*h-alpha*a*(b+c)-alpha*a*(d+e);
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+1)=-alpha*b*d-alpha*b*e-alpha*b^2-alpha*a*b-alpha*b*c-alpha*b*h;
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+2)=-alpha*c*e-alpha*c*d-alpha*c^2-alpha*a*c-alpha*b*c-alpha*c*h;
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+3)=0;
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+4)=(1/2)*alpha*m*Pr*(b*e+c*d)+alpha*(1-(1-m/2)*Pr)*a*(d+e)+alpha*a*h+(1-m/2)*alpha*Pr*(b*h+c*h)+alpha*a*(b+c)+alpha*a^2+alpha*m*Pr*b*c;
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+5)=(1-m/2)*alpha*Pr*a*d+alpha*(1-Pr)*b*e+alpha*b*d+alpha*b^2+alpha*a*b+alpha*(1-m*Pr)*b*c+alpha*(1-(1-m/2)*Pr)*b*h+(1-m)*alpha*Pr*c*d;
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+6)=(1-m/2)*alpha*Pr*a*e+alpha*(1-Pr)*c*d+alpha*c*e+alpha*c^2+alpha*a*c+alpha*(1-m*Pr)*b*c+alpha*(1-(1-m/2)*Pr)*c*h+(1-m)*alpha*Pr*b*e;
    dx(idx_genotypes_in_meiosis(i)+7*(idx_genotypes_in_meiosis(i)-1)+7)=(1/2)*alpha*m*Pr*(b*e+c*d)+alpha*m*Pr*b*c;

    end
end
    

 for i=1:length(idx_genotypes_not_in_meiosis)        
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)) = r*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1));  
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1) = r*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1);  
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2) = r*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2); 
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+3) = 0;           
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+4) = 0;    
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+5) = 0;  
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+6) = 0;  
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+7) = 0;
       
       for j=1:8*S       
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)) = dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1))-(r/k)*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1))*x(j);  
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1) = dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1)-(r/k)*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1)*x(j);
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2) = dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2)-(r/k)*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2)*x(j);
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+3) = 0; 
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+4) = 0;
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+5) = 0; 
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+6) = 0; 
       dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+7) = 0;
       end
       
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1))=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1))-2*muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1));
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1)+muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1))-muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1);
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2)+muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1))-muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2);
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+3)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+3)+muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+1)+muH*x(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+2);         
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+4)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+4);         
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+5)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+5);
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+6)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+6);
         dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+7)=dx(idx_genotypes_not_in_meiosis(i)+7*(idx_genotypes_not_in_meiosis(i)-1)+7);
 end