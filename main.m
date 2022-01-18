%%%%%%%%-----------------------%%%%%%%%%%%
clc
clear all
close all
%%%%%%%%-----------------------%%%%%%%%%%%


%%%%%%%%-----------------------%%%%%%%%%%%
d=2; %% Dimension 2 or 3 only %%
sf=1.1; %%Scale factor%%
gausianparameter=power(10,-5); %%Gaussian Parameter%%
nsample=20; %%number of points
discre=power(10,-10); %%discrepancy
p1=eye(d); %% Initial Point

%%%%%%%%------------Setup-----------%%%%%%%%%%%
lam=0.1;
interaciones=20;
alamar=0;
s=size(p1);
dimen=s(2);
p1=funcionnorma(p1,dimen);
p=reshape(p1,1,[]);
s=size(p);
conta=2;
counter=0;
fante=0;
volumenesfe=power(pi,dimen/2)/(gamma(1+(dimen/2)));
maxpoint=p1;

%%%%%%%%-----------------------%%%%%%%%%%%


[x1, x2]=lenum(p1);
eva=norm(x1);
valormax=power((eva/2),dimen)*volumenesfe;
format long
nsolsp=pi/(2*sqrt(3));
ce={};
gua={};

while alamar==0;
 
    
    p=reshape(p1,1,[]);

    sigma=gausianparameter*gausianparameter*eye(s(2));
    R = mvnrnd(p,sigma,nsample);
    inte=0*p1;
    
    %fx= funcionopti(p1,dimen);
    [x1, x2]=lenum(p1);
    fx=norm(x1);
    fx=power((fx/2),dimen)*volumenesfe;
    


    
  for j=1:1:nsample
     
      
      ytest=reshape(R(j,:),dimen,[]);
      
   ytest=funcionnorma(ytest,dimen);
   [y1, y2]=lenum(ytest);
    fy=norm(y1);
    fy=power((fy/2),dimen)*volumenesfe;
   
   
      
      
      %fy= funcionopti(ytest,dimen);

 
  
      inte=inte+(ytest-p1)*(fy-fx);
  
      
      
  end
  
  inte=inte/(nsample*gausianparameter*gausianparameter);
  
  p1=p1+lam*inte;
 
  p1=funcionnorma(p1,dimen);
 

  %eva=funcionopti(p1,dimen);
  [x1, x2]=lenum(p1);
    eva=norm(x1);
  estimativa=power((eva/2),dimen)*volumenesfe;





if( valormax < estimativa)
     maxpoint=p1;
    valormax=estimativa;

end

if( interaciones < counter )
counter=0;
p1=maxpoint;
lam=lam/sf;
end

 if( abs(estimativa-fante) <discre)
  alamar=1;
 end
 
 
 ce{conta-1}=abs(nsolsp-estimativa);
  gua{conta-1}=p1;
 conta=conta+1;
 counter=counter+1;
  
fante=estimativa;
end
  
  
%%%%------------plot-------------%%%%

if(d==2)
    dibujar2(p1);
    
end
 
 
 if(d==3)
    dibujar(p1);
    
 end
 
%%%%------------Error-------------%%%%

cek = Error(ce);
