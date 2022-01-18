function [solum] = funcionnorma(M,dimen)
    
 deter=det(M);
 
%aba=abs(deter);

vec=sqrt(sum(M.^2,2));

if (deter < 1)
    [vecmin,I] = min(vec);
    M(:,I)=M(:,I)./deter;
    
end
    
if (deter > 1)
    [vecmin,I] = max(vec);
    M(:,I)=M(:,I)./deter;
    
end
    

 solum=M;
 
end