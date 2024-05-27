function dx = ODE_sex_arbitrary(t,x,r0,r1,r2,muM,alpha,k,idx_genotypes_in_meiosis,idx_genotypes_not_in_meiosis,S,Pr)


dx=zeros(4*S,1);


if ~isempty(idx_genotypes_in_meiosis)

    
    for i=1:length(idx_genotypes_in_meiosis)
    d=0; e=0; h=0;    

    a=x(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1));
    b=x(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1)+1);
    c=x(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1)+2);
    

      for j=1:length(idx_genotypes_in_meiosis)

         if j~=i
         xx=x(idx_genotypes_in_meiosis(j)+3*(idx_genotypes_in_meiosis(j)-1)+1);
         d=d+xx;
         yy=x(idx_genotypes_in_meiosis(j)+3*(idx_genotypes_in_meiosis(j)-1)+2);
         e=e+yy;    
         zz=x(idx_genotypes_in_meiosis(j)+3*(idx_genotypes_in_meiosis(j)-1));
         h=h+zz;
         end

      end
    dx(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1))=(1/3)*alpha*Pr*(b*e+c*d)+(2/3)*alpha*Pr*b*c+(2/3)*alpha*Pr*(b+c)*h-(2/3)*alpha*Pr*a*(d+e);
    dx(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1)+1)=(2/3)*alpha*Pr*a*d-(2/3)*alpha*Pr*b*h-alpha*Pr*b*e+(1/3)*alpha*Pr*c*d-(2/3)*alpha*Pr*b*c;
    dx(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1)+2)=(2/3)*alpha*Pr*a*e-(2/3)*alpha*Pr*c*h-alpha*Pr*c*d+(1/3)*alpha*Pr*b*e-(2/3)*alpha*Pr*b*c;
    dx(idx_genotypes_in_meiosis(i)+3*(idx_genotypes_in_meiosis(i)-1)+3)=(1/3)*alpha*Pr*(b*e+c*d)+(2/3)*alpha*Pr*b*c;

    end

end
    

 for i=1:length(idx_genotypes_not_in_meiosis)        
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)) = r0*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1));
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1) = r1*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1);
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2) = r1*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2);
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3) = r2*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3);
       
       for j=1:4*S       
       dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)) = dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1))-(r0/k)*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1))*x(j);
       dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1) = dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1)-(r1/k)*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1)*x(j);
       dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2) = dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2)-(r1/k)*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2)*x(j);
       dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3) = dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3)-(r2/k)*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3)*x(j);
       end
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1))=dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1))-2*muM*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1));
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1)=dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1)+muM*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1))-muM*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1);
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2)=dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2)+muM*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1))-muM*x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2);
         dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3)=dx(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+3)+muM*(x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+1)+x(idx_genotypes_not_in_meiosis(i)+3*(idx_genotypes_not_in_meiosis(i)-1)+2));

 end