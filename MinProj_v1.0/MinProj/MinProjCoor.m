%%
% MinProjCoor.m
% Computes minimal projections onto coordinate subspaces
% by solving the linear optimization program
% outlined in the article
% COMPUTATION OF MINIMAL PROJECTIONS AND EXTENSIONS
% by S. Foucart

% Find solutions to
% minimize ||P||_{1 -> 1} or ||P}|_{infty -> infty}
% over all projections P from R^n onto a given m-dimensional subspace
%
% Usage: [projCst,minProj,status] = MinProjCoor(U,q,...)
%
% U: n-by-m matrix whose columns span the m-dimensional subspace
% q: the index of the norm, either 1 or inf
% optional inputs: parameters controlling CVX
%
% projCst: the projection constant of the subspace of R^n relative to the q-norm
% minProj: n-by-n matrix (relative to the canonical basis of R^n) of a minimal projection
% status: reliability of the result (duplicates cvx_status)
%
% Written by Simon Foucart in May 2014
% Last updated in December 2014
% Send comments to simon.foucart@centraliens.net

function [projCst, minProj, status] = MinProjCoor(U, q, varargin)

% parameters controlling CVX
% set cvx_quiet 
loc = find(strcmpi(varargin,'quiet'));
if any(loc)
  cvx_quiet(varargin{loc+1});
else
  cvx_quiet(true);
end
% set cvx_precision ('low', 'medium', 'default', 'high', 'best')
loc = find(strcmpi(varargin,'precision'));
if any(loc)
  cvx_precision(varargin{loc+1});
else
  cvx_precision('default');
end
% choose the sdp solver ('SDPT3', 'SeDuMi', or 'Mosek')
loc = find(strcmpi(varargin,'solver'));
if any(loc)
  cvx_solver(varargin{loc+1});
end

% definition of a matrix associated with the orthogonal complement of U
[n,m] = size(U);
[Q,~] = qr(U);
Utilde=Q(:,m+1:n);

cvx_begin

% introduction of the optimization variables
variable d;
variable vecP(n*n);
variable c(n*n);

% minimization of the objective function
minimize d

% formulation of the "projection" constraints
subject to
kron(U',eye(n))*vecP == U(:);
kron(Utilde',Utilde')*vecP == 0;

% formulation of the slack constraints
c+vecP >= 0;
c-vecP >= 0;
if q==1
  for j=1:n
    sum( c((j-1)*n+1:j*n) ) <= d;
  end
end
if q==inf
  for i=1:n
    sum( c(i+n*(0:n-1)) ) <= d;
  end
end

cvx_end

% return the outputs
projCst = cvx_optval;
minProj = reshape(vecP,n,n);
status = cvx_status;
if ~(strcmp(status,'Solved'))
  warning(strcat('the optimization status is', 32, status))
end

end