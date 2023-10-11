function [rs ks]=genfig8(rin,rout,height,Nloop);
[zval,wtz]=lgwt(2,0,height);
[rval,wtr]=lgwt(20,0,rout);
wtr=wtr.*rval;
phi=2*pi/36:2*pi/36:2*pi;
wtphi=ones(size(phi));
[rval phi zval]=ndgrid(rval(:),phi(:),zval(:));
[wtr wtphi wtz]=ndgrid(wtr(:),wtphi(:),wtz(:));
wt=wtr(:).*wtphi(:).*wtz(:)/height;
wt(rval(:)<rin)=Nloop.*wt(rval(:)<rin);
wt(rval(:)>=rin)=(Nloop-Nloop/(rout-rin)*(rval(rval(:)>=rin)-rin)).*wt(rval(:)>=rin);
rs(:,1)=rval(:).*cos(phi(:));
rs(:,2)=rval(:).*sin(phi(:));
rs(:,3)=zval(:);
ks=zeros(size(rs));
ks(:,3)=wt(:);
rs(:,1)=rs(:,1)-rout;
rs(end+1:2*end,:)=rs(:,:);
rs(end/2+1:end,1)=rs(end/2+1:end,1)+2*rout;
ks(end+1:2*end,:)=-ks(:,:);
% quiver3(rs(:,1),rs(:,2),rs(:,3),ks(:,1),ks(:,2),ks(:,3));
% hold on
end
function [x,w]=lgwt(N,a,b)

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
