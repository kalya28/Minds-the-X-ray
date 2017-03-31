function [h1,h]=c_means(im)
[x y z]= size(im);

if( z==3)
    im=rgb2gray(im);
end

a=im;
b=im2double(a);

[m n]=size(a);

k=2;
cc1=0.9;
cc2=0.3;
c=cat(3,b,b);
iteration=0;


while(iteration<2)
    iteration=iteration+1
    
    d=repmat(cc1,m,n);
    e=repmat(cc2,m,n);
    if iteration==1 
        test1=d; test2=e;
    end
    f=cat(3,d,e);
    
    ree=repmat(0.000001,m,n); 
    ree1=cat(3,ree,ree);
    
    distance=c-f;
    distance=distance.*distance+ree1;
    
    daoShu=1./distance;
    
    daoShu2=daoShu(:,:,1)+daoShu(:,:,2);
    distance1=distance(:,:,1).*daoShu2;
    u1=1./distance1;
    distance2=distance(:,:,2).*daoShu2;
    u2=1./distance2;
      
    ccc1=sum(sum(u1.*u1.*b))/sum(sum(u1.*u1));
    ccc2=sum(sum(u2.*u2.*b))/sum(sum(u2.*u2));
   
    tmpMatrix=[abs(cc1-ccc1)/cc1,abs(cc2-ccc2)/cc2];
    pp=cat(3,u1,u2);
    
    for i=1:m
        for j=1:n
            if max(pp(i,j,:))==u1(i,j)
                g(i,j)=1;
           
            else
                g(i,j)=2;
            end
   if max(tmpMatrix)<0.00001
         break;
  else
         cc1=ccc1;
         cc2=ccc2;
        
  end

 for i=1:m
       for j=1:n
            if g(i,j)==2
            h(i,j)=125;
                 else
            h(i,j)=12;
       end
    end
end

end

for i=1:m
    for j=1:n
         if g(i,j)==2
            h(i,j)=185;
             else
             h(i,j)=10;
    end
  end
end 




