function [solum] = dibujar2(M)
    
[a,b]=lenum(M);
r=norm(a)/2;


figure

for n= -3:1:3
for m= -3:1:3
  
    v=[ n ; m ];
    p=transpose(M*v);

    p1=p(1);
    p2=p(2);

    

ang=0:0.01:2*pi; 
xp=p1+r*cos(ang);
yp=p2+r*sin(ang);

plot(xp,yp,'k');
hold on


end
end

axis equal




end