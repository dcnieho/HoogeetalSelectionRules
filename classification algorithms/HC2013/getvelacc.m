function [V, A, C, tim] = getvelacc(Y,dn,freq)
% GETVELACC(Y, dn, freq) ;
%
% Computes the instantaneous velocity and acceleration of position samples Y using a
% polynomial fit on 2*dn+1 points.
%
% [V, A] = GETVEL returns velocity V and acceleration A
% [V, A, C] = GETVEL also returns fitted position on Y
% [V,A,C,tim] = GETVEL returns time consumed by program in seconds
%
% GETVEL(Y) uses dn=3 and freq = 250 ; (EyeLink Sample Freq = 250 Hz).
%
%
% Y(j) = C(i) + V(i) * X(j) + A(i) * X(j)^2
% in which X(j) is position at sample i-dn:i+dn
% V(i) and A(i) are the fitted velocity and acceleration at point i
%
% 2000 JN van der Geest

% Theoretically, we solve:
% SUM(y-y') minimal, in which y' = a + bx+cx^2
% gives three equations:
%   [  Na         - SUM(y)     + b*SUM(x)   + c*SUM(x^2) = 0
%   [  b*SUM(x^2) - SUM(xy)    + a*SUM(x)   - c*SUM(x^3) = 0
%   [  c*SUM(x^4) - SUM(y*x^2) + a*SUM(x^2) + b*SUM(x^3) = 0
%
% in which x=-dn:dn, and therefore N = 2*dn+1, SUM(x) = SUM(x^3) = 0,
% yields:
% b = SUM(x*y) / SUM(x^2)                                             = 'Velocity'
% c = [N*SUM(y*x^2) - SUM(x^2)*SUM(y)] / [N*SUM(x^4) - (SUM(x^2))^2]  = 'Acceleration'
% a = [SUM(y) - c * SUM(x^2)] / SUM(x^4)                              = 'Position'

%----------------------------------------
% size to reshape V,A and C
if nargin==1,
   dn = 3 ;
   freq = 250 ;
end
if nargin==2,
   freq = 250 ;
end

tim = cputime ;
[nry ncy] = size(Y) ;

% make a row vector
Y = Y(:)' ;
N = length(Y) ;
x = -dn:dn ; x = x(:) ; Nx = length(x) ; % = 2*dn + 1

% create an index matrix
% x2(:,j) = x(j) ;
x2 = repmat(x,1,N) ; 
%ind(i,j) = i+x(j) ;
ind = repmat(1:N,Nx,1) + x2  ;

% create a data matrix, index falling outside range of Y will be nans
q = (ind<1|ind>N) ; 
ind(q) = 1 ;
% Y(i,j) = Y(i+x(j)) 
Y = Y(ind) ;
Y(q) = nan ;

% now we treat x2 and y columnwise and fit 2nd degree polynomial
% S.. vectors are row vectors with the same length as Y
x2 = x2 / freq ;
Sx2 = sum(x2 .^ 2) ;
Sx4 = sum(x2 .^ 4) ;
Sy = sum(Y) ;
Sx2y = sum(x2 .* x2 .* Y) ;
Sxy = sum (x2 .* Y) ;

% calculate polynomial coefficients element by element
% and return in the same shape as Y ;
V = reshape((Sxy ./ Sx2),nry,ncy) ;
if nargout > 1,
   A = reshape(((Nx .* Sx2y) - (Sx2 .* Sy)) ./ ((Nx .* Sx4) - (Sx2 .* Sx2)),nry,ncy) ;
   if nargout>2,
      C = reshape(((Sy - (A(:)' .* Sx2)) ./ Nx),nry,ncy)  ;
   end
end
tim = cputime - tim ;

