function [I_fin] = VPI_dicom(I_in,theta,s)
% prodotto tensoriale 3d
% s scalare fattore di scala
[n1,n2,n3]=size(I_in);
a=length(s); special=0;
switch a
    case 1 %s= real scale factor
      %N1=fix(s*n1); N2=fix(s*n2);
      N1=round(s*n1); N2=round(s*n2);N3=round(s*n3);
      if s<1 && rem(1/s,2)==1
          special=1; ss=1/s;
      end
   case 3 %s=[Nrow, Ncol];
      N1=s(1);N2=s(2);N3=s(3);
      if N1<n1 && rem(n1/N1,2)==1 && N2<n2 && rem(n2/N2,2)==1 && N3<n3 && rem(n3/N3,2)==1
          special=2; ss=[n1/N1, n2/N2, n3/N3];
      end
    otherwise
      special=3; %error indicator about the last input parameter
end
if theta<0 || theta>1
    special=4; %error indicator about the second input parameter
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch special
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Special cases: Downscaling with odd scale factors or input errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1 
      I_fin=I_in((ss*(2*(1:N1)-1)+1)/2,(ss*(2*(1:N2)-1)+1)/2,:);
  case 2
      I_fin=I_in((ss(1)*(2*(1:N1)-1)+1)/2,(ss(2)*(2*(1:N2)-1)+1)/2,:);
  case 3
      fprintf(' some parameter not adequate\n'); 
end




%% Computing the new size N1xN2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % N1=round(s*n1); N2=round(s*n2);N3=round(s*n3);
    eta=(2*(1: N1)-1)*pi/(2*N1); 
    csi=(2*(1: N2)-1)*pi/(2*N2);  
    tau=(2*(1: N3)-1)*pi/(2*N3);
%--------------------------------
% Computing the basis polynomials
%---------------------------------
  if theta==0 % LCI method
    n1MAX=n1-1; n2MAX=n2-1; n3MAX=n3-1;
    T1=cos((0:n1MAX)'.*eta)*sqrt(2/n1);T1(1,:)=sqrt(1/n1);
    T2=cos((0:n2MAX)'.*csi)*sqrt(2/n2);T2(1,:)=sqrt(1/n2);
    T3=cos((0:n3MAX)'.*tau)*sqrt(2/n3);T3(1,:)=sqrt(1/n3);
    lx=idct(T1); ly=idct(T2); lz=idct(T3);
  else
    m1=fix(theta*n1); m2=fix(theta*n2);m3=fix(theta*n3);
    n1MAX=n1+m1-1; n2MAX=n2+m2-1;n3MAX=n3+m3-1;
    q1=zeros(n1,N1); q2=zeros(n2,N2);q3=zeros(n3,N3);
    T1=cos((0:n1MAX)'.*eta)*sqrt(2/n1);T1(1,:)=sqrt(1/n1);
    T2=cos((0:n2MAX)'.*csi)*sqrt(2/n2);T2(1,:)=sqrt(1/n2);
    T3=cos((0:n3MAX)'.*tau)*sqrt(2/n3);T3(1,:)=sqrt(1/n3);
    n1meno=n1-m1+1;n2meno=n2-m2+1; n3meno=n3-m3+1;
    q1(1:n1meno,:)=T1(1:n1meno,:);q2(1:n2meno,:)=T2(1:n2meno,:);
    q3(1:n3meno,:)=T3(1:n3meno,:);
    z=[m1-1:-1:1];
    q1(n1meno+1:n1,:)=((m1+z')/(2*m1)).*T1(n1-z+1,:)+...
        ((z'-m1)/(2*m1)).*T1(n1+z+1,:);
    z=[m2-1:-1:1];
    q2(n2meno+1:n2,:)=((m2+z')/(2*m2)).*T2(n2-z+1,:)+...
        ((z'-m2)/(2*m2)).*T2(n2+z+1,:);
    z=[m3-1:-1:1];
    q3(n3meno+1:n3,:)=((m3+z')/(2*m3)).*T3(n3-z+1,:)+...
        ((z'-m3)/(2*m3)).*T3(n3+z+1,:);
    lx=idct(q1); ly=idct(q2);lz=idct(q3); 
  end
  % in ogni caso 
  lz=lz';% dimensioni N3xn3
  %-----------------------------
% Computing the resized image
%----------------------------

%size(lz)
    I_fin=zeros(N1,N2,N3);
    I_in=double(I_in);% 
    
    for j=1:n3
        mat(:,:,j)=(lx'*I_in(:,:,j))*ly;
    end

    for i=1:N3
        for j=1:n3
           I_fin(:,:,i)=I_fin(:,:,i)+mat(:,:,j)*lz(i,j);
        end
    end
    %I_fin=uint8(I_fin);
end