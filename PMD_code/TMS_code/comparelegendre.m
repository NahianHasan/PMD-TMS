function [outP,outL,outExact]=comparelegendre(funct,evalpts,ord)
outExact=funct(evalpts);
outP=zeros(size(evalpts));
outL=outP;
[x,w]=lgwt(ord+1,-1,1);
%make legendre interpolation
legendv=zeros([ord+1 ord+1]);
for i=1:ord+1
    legendv(:,i)=legendre_polynomial (ord,x(i));
end
samplepts=funct(x);
inner_prod=2./(2*(0:ord)+1);
w=(w(:).*samplepts(:));
coeff=zeros([ord+1 1]);
for i=1:ord+1
coeff(i)=inner_prod(i)^-1*(legendv(i,:)*w);
end
for i=1:numel(evalpts)
    legends=legendre_polynomial (ord,evalpts(i));
outP(i)=coeff(:)'*legends(:);
end

%%compute lagrange value
    outL=lagrange(reshape(evalpts,[1,numel(evalpts)]),reshape(x,[1,ord+1]),reshape(samplepts,[1,ord+1]))';
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


function y=lagrange(x,pointx,pointy)
%
%LAGRANGE   approx a point-defined function using the Lagrange polynomial interpolation
%
%      LAGRANGE(X,POINTX,POINTY) approx the function definited by the points:
%      P1=(POINTX(1),POINTY(1)), P2=(POINTX(2),POINTY(2)), ..., PN(POINTX(N),POINTY(N))
%      and calculate it in each elements of X
%
%      If POINTX and POINTY have different number of elements the function will return the NaN value
%
%      function wrote by: Calzino
%      7-oct-2001
%
n=size(pointx,2);
L=ones(n,size(x,2));
if (size(pointx,2)~=size(pointy,2))
   fprintf(1,'\nERROR!\nPOINTX and POINTY must have the same number of elements\n');
   y=NaN;
else
   for i=1:n
      for j=1:n
         if (i~=j)
            L(i,:)=L(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
         end
      end
   end
   y=0;
   for i=1:n
      y=y+pointy(i)*L(i,:);
   end
end
end
function l = legendre_polynomial ( p, x )

%*****************************************************************************80
%
%% LEGENDRE_POLYNOMIAL evaluates the Legendre polynomials L(0:P)(X).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 May 2008
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Milton Abramowitz and Irene Stegun,
%    Handbook of Mathematical Functions,
%    US Department of Commerce, 1964.
%
%    Daniel Zwillinger, editor,
%    CRC Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996.
%
%  Parameters:
%
%    Input, integer P, the highest evaluation degree.
%
%    Input, real X, the evaluation point.
%
%    Output, real L(1:P+1), the Legendre polynomials of order 0 through P at X.
%
  if ( p < 0 )
    l = [];
    return
  end
%
%  Allocate space.
%
  l(1:p+1) = 0.0;
%
%  Apply recursion.
%
  l(1) = 1.0;

  if ( 1 <= p )

    l(2) = x;
 
    for i = 2 : p
 
      l(i+1) = ( ( 2 * i - 1 ) * x * l(i)   ...
               - (     i - 1 ) *     l(i-1) ) ...
               / (     i     );
 
    end

  end
%
%  Normalize.
%

  return
end