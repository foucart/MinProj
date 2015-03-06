%%
% BlatterCheney.m
% Computes the projection constant of a hyperplane according the formula in
% MINIMAL PROJECTIONS ON HYPERPLANES IN SEQUENCE SPACES
% by Blatter and Cheney.

% Find solutions to
% minimize ||P||_{1 -> 1} or ||P||_{infty -> infty\}
% over all projections P from R^n onto a given hyperplane
%
% Usage: projCst = BlatterCheney(f,q)
%
% f: a vector in R^n orthogonal to the hyperplane under consideration
% q: the index of the norm, either 1 or Inf
%
% projCst: the projection constant of the hyperplane relative to the q-norm
%
% Written by Simon Foucart in June 2014
% Last updated in December 2014
% Send comments to simon.foucart@centraliens.net

function projCst = BlatterCheney(f,q)

if q == Inf
   f = f/norm(f,1);
   if norm(f,Inf) <1/2
      alpha = abs(f) ./ (1-2*abs(f));
      projCst = 1 + 1./sum(alpha);
   else
      projCst = 1 ;
   end
end

if q == 1
   f = sort(abs(f),'descend');
   f = f/f(1);
   if f(3) == 0
      projCst = 1;
   end
   if f(3) > 0
      a = cumsum(f);
      b = cumsum(1./f);
      c = min(f(2:end).*b(1:end-1),a(1:end-1));
      v = find( c - (1:length(c))' + 2 < 0);
      if isempty(v)
         disp('Error: case not covered')
      else
         n = min(v);
         projCst = 1+2/...
            ( a(n)*b(n)/(n-2) - n + (b(n)/(n-2)-1/f(n))*max(0,n-2-a(n)) );
      end
   end
   
end

end