%%
% ExMinExt.m
% Computes the extension constant of a specific linear map 
% by solving the linear or semidefinite optimization programs
% as outlined in the last section of the article
% COMPUTATION OF MINIMAL PROJECTIONS AND EXTENSIONS
% by S. Foucart

% Find solutions to
% minimize |||P|||
% over all extensions from \ell_infty of a map A onot P_\infty^l 
%
% Usage: [extCst,minExt] = ExMinExt(l,pts)
%
% l: the dimension of the polynomial space P_\infty^l
% pts: a column vector containing the points t_1<...t_n defining A
%
% extCst: the extension constant, i.e., teh value of the minimium
% minExt: a minimal extension, i.e., a minimizer
%
% Written by Simon Foucart in December 2014
% Send comments to simon.foucart@centraliens.net

function [extCst,minExt] = ExMinExt(l,pts)

n = length(pts);
U = zeros(n,l);
for j = 1:l
  Tj = chebpoly(j-1);
  U(:,j) = Tj(pts);
end
Il = eye(l);

cvx_quiet true;

cvx_begin sdp

variable d;
variable P(l,n);
variable Q(l,l,2^n) symmetric semidefinite;

minimize d

subject to
kron(U',Il)*P(:) == Il(:);
for h = 1:2^n
  eps_aux = dec2bin(2^n+h-1)-'0';
  eps = 2*(eps_aux(2:end)')-1;   %this creates the h-th sequence of -1 and 1
  % k>=1
  for k = 1:l-1
   P(k+1,:)*eps/2 == sum( diag(Q(:,:,h),-k) );
  end
  % k=0
  d+P(1,:)*eps == sum( diag(Q(:,:,h)) );
end

cvx_end

extCst = cvx_optval;
minExt =  P;

end