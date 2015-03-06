%% 
% MinProjPoly.m
% Computes minimal projections onto polynomial subspaces
% by solving the linear or semidefinite optimization programs
% outlined in the article
% COMPUTATION OF MINIMAL PROJECTIONS AND EXTENSIONS
% by S. Foucart

% Find solutions to
% minimize ||P||_{N -> N} or ||P||_{infty -> N}
% over all projections P from P^n onto a given m-dimensional subspace
% and return lower and upper bounds for the minimum of ||P||_{infty -> infty}
%
% Usage: [bounds,proj,status] = MinProjPoly(basis,n,nu,...)
%
% basis: a cell whose elements are polynomials (as chebfuns) spanning the m-dimensional subspace
% n: the dimension of the polynomial superspace P^n
% nu : the discretization index
% optional inputs: 'LP' or 'SDP' for the optimization program being solved (default: 'SDP') 
% and parameters controlling CVX
%
% bounds: lower and upper bounds for the projection constant of the subspace in P^n
% proj: n-by-n matrix of a pseudo minimal projection
% (relative to the Lagrange basis for 'LP' and to the Chebyshev basis of 'SDP')
% status: reliability of the result (duplicates cvx_status)
%
% Written by Simon Foucart in December 2014
% Send comments to simon.foucart@centraliens.net

function [bounds, proj, status] = MinProjPoly(basis, n, nu, varargin)

% choose the approximation problem to be solved, i.e.,
% linear program (LP) or semidefinite program (SDP)
loc = find(strcmpi(varargin,'LP'));
if any(loc)
  LP = true;
  SDP = false;
end
loc = find(strcmpi(varargin,'SDP'));
if any(loc)
  LP = false;
  SDP = true;
end
if ( ~exist('LP') && ~exist('SDP'))
  LP = false;
  SDP = true;
end

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

% start preparing data for both the LP and SDP cases
m = length(basis);
if max(cellfun(@length,basis)) > n
  error('subspace spanned by first inputs not contained in space with second input as dimension')
end
N = (n-1)*nu + 1;
tN = chebpts(N);


%% The linear program

if LP
  
  % continue preparing the data for LP
  U = zeros(n,m);
  tn = chebpts(n);
  for j = 1:m
    uj = basis{j};
    U(:,j) = uj(tn);
  end
  [Q,~] = qr(U);
  Utilde = Q(:,m+1:n);
  L = chebfun.lagrange(tn);
  B = zeros(N,n);
  for i = 1:N
    for j = 1:n
      Lj = L(:,j);
      B(i,j) = Lj(tN(i));
    end
  end
  Insert = zeros(N,n);
  In = eye(n);
  for i = 1:n
    Insert((i-1)*nu+1,:) = In(i,:);
  end
  
  cvx_begin
  % introduction of the optimization variables
  variable d;
  variable P(n,n);
  variable X(N,N);
  variable Y(N,N);
  variable Z(N,N);
  
  % minimization of the objective function
  minimize d
  
  % formulation of the "projection" constraints
  subject to
  kron(U',eye(n))*P(:) == U(:);
  kron(Utilde',Utilde')*P(:) == 0;
  
  % formulation of the other constraints
  X >= 0;
  Y >= 0;
  X - Y + Z - Insert*B'*Z == -Insert*P'*B';
  for i =1:N
    sum(X(:,i)) + sum(Y(:,i)) <= d;
  end
  
  cvx_end
  
  % return the output
  bounds = zeros(1,2);
  bounds(1) = cvx_optval*(1-pi/2/nu);
  bounds(2) = cvx_optval/(1-pi/2/nu);
  proj = P; 
  status = cvx_status;
  if ~(strcmp(status,'Solved'))
    warning(strcat('the optimization status is', 32, status))
  end
  
end

%% The semidefinite program

if SDP
  
  % continue preparing the data for SDP
  U = zeros(n,m);
  for j = 1:m
    uj = basis{j};
    coeffsj = chebcoeffs(uj);
    U(1:length(coeffsj),j) = coeffsj;
  end
  [Q,~] = qr(U);
  Utilde = Q(:,m+1:n);
  V = zeros(N,n);
  for j = 1:n
    Tj = chebpoly(j-1);
    V(:,j) = Tj(tN);
  end
  
  cvx_begin sdp
  % introduction of the optimization variables
  variable d;
  variable P(n,n);
  variable Z(n,n,N) symmetric toeplitz;
  
  % minimization of the objective function
  minimize d
  
  % formulation of the "projection" constraints
  subject to
  kron(U',eye(n))*P(:) == U(:);
  kron(Utilde',Utilde')*P(:) == 0;
  
  % formulation of the other constraints
  for i = 1:N
    Z(:,:,i) >= +toeplitz(V(i,:)*P);
    Z(:,:,i) >= -toeplitz(V(i,:)*P);
    Z(1,1,i) <= d;
  end
  
  cvx_end
  
  % return the output
  bounds = zeros(1,2);
  bounds(1) = cvx_optval;
  bounds(2) = cvx_optval/(1-pi/2/nu);
  proj = P;
  status = cvx_status;
  if ~(strcmp(status,'Solved'))
    warning(strcat('the optimization status is', 32, status))
  end
 
end

end