function [solum] = dibujar(M)
    

[x,y,z] = sphere;
p1=M(:,1);
p2=M(:,1);
p3=M(:,1);


[a,b]=lenum(M);
r=norm(a)/2;
x=r*x;
y=r*y;
z=r*z;
figure

for n= -2:1:2
for m= -2:1:2
for k= -2:1:2
  
    v=[ n ; m ; k];
    p=M*v;
    p1=p(1);
p2=p(2);
p3=p(3);

surf(x+p1,y+p2,z+p3)
colormap('summer');
hold on


end
end
end


axis equal

end