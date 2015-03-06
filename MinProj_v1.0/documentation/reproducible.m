%% Reproducible file accompanying the manuscript
% Computattion of minimal projections and extensions
% by Simon Foucart

%%
% This file generates the reproducible experiments mentioned in [1]
% that illustrate how to use the MATLAB package |MinProj| (version 1.0).
% The package relies heavily on CVX [2] and Chebfun [3],
% which are included in the downloadable folder.
% CVX and Chebfun need to be activated by typing |MinProj_setup|
% at MATLAB's prompt from the folder MinProj_v1.0.


%% Verification of exact formulas for coordinate spaces

%%
% The following code generates a number of random hyperplanes
% and compares the values of the projection constants relative to $\ell_\infty^n$
% predicted by Blatter-Cheney's formula [4]
% and computed via the function MinProjCoor
% (hence it validates the implementation of MinProjCoor).

n = 10;          % dimension of the space V=\ell_\infty^n
nTests = 50;     % number of random subspaces generated
equalities = zeros(nTests,1);
tic;
for k = 1:nTests;
  f = randn(n,1);
  p_theor = BlatterCheney(f,Inf);      % theoretical projection constant of f^\perp
  [Q,R] = qr(f); U = Q(:,2:n);         % U is the space f^\perp
  p_comp = MinProjCoor(U,Inf);         % computed projection constant of f^\perp
  equalities(k) = ( norm(p_theor-p_comp)<1e-5 );   % record 1 for equality
end
t = toc;
%display the result of the test
if sum(equalities) == nTests
  disp(strcat('Theoretical and computed values always agreed (',...
    num2str(round(t)), 'sec)'))
end

%%
% The same verification is now performed for the projection constant
% relative to $\ell_1^n$ rather than to $\ell_\infty^n$.

n = 10;          % dimension of the space V=\ell_\infty^n
nTests = 50;     % number of random subspaces generated
equalities = zeros(nTests,1);
tic;
for k = 1:nTests;
  f = randn(n,1);
  p_theor = BlatterCheney(f,1);       % theoretical projection constant of f^\perp
  [Q,R] = qr(f); U = Q(:,2:n);        % U is the space f^\perp
  p_comp = MinProjCoor(U,1);          % computed projection constant of f^\perp
  equalities(k) = ( norm(p_theor-p_comp)<1e-5 );   % record 1 for equality
end
t = toc;
%display the result of the test
if sum(equalities) == nTests
  disp(strcat('Theoretical and computed values always agreed (',...
    num2str(round(t)), 'sec)'))
end

%%
% The following verification concerns Lewicki's formula
% for 2-codimensional subspaces of $\ell_\infty^n$ [5].
% Note that not all randomly generated subspaces will satify the conditions
% of applicability of the formula, so a test of these conditions is included.

nTests = 5000;     % number of random 2-dimensioanl subspaces generated
nApplicable = 0;
equalities = [];
tic;
for k = 1:nTests;
  f = [rand(1,1); 0 ; rand(2,1)]; f = f/norm(f,1);
  g = [0; rand(3,1)]; g = g/norm(g,1);
  [p_theor, isApplicable] = Lewicki(f,g,0);   %theoretical projection constant of [f,g]^\perp
  if isApplicable == 1
    nApplicable = nApplicable + 1;
    [Q,R] = qr([f g]); U = Q(:,3:4);   %U is the space [f,g]^\perp
    p_comp = MinProjCoor(U,Inf);       %computed projection constant of a^\perp
    equalities = [equalities, ( norm(p_theor-p_comp)<1e-5 )];
  end
end
t = toc;
%display the result of the test
disp(strcat('Out of', 32, num2str(length(equalities)), ...
  ' applicable cases, theoretical and computed values aggreed ', 32, ...
  num2str(sum(equalities)), ' times (', num2str(round(t)), 'sec)'))


%% Exploration of expected projection constants in $\ell_\infty^n$

%%
% THIS CELL IS NOT INTENTED TO BE EVALUATED! IT WOULD TAKE ABOUT 20 DAYS!
% The following code generates random $m$-dimensional subspaces of $\ell_\infty^n$
% and computes the average value of the projection constants
% for different values of $m$ and $n$.

try load('AvProj_default.mat')
catch
  m_max = 16;
  n_max = 5/2*m_max;
  NbTests = 100;
  Proj_cst = zeros(m_max,n_max);
  tic;
  for m = 1:m_max
    for n = m+1:n_max
      for k = 1:NbTests
        U = randn(n,m);
        Proj_cst(m,n) = Proj_cst(m,n) + MinProjCoor(U,Inf);
      end
    end
  end
  Proj_cst = Proj_cst/NbTests;
  t = toc;
  save('AvProj.mat')
end

%%
% Here is a graphical visualization of the results

try load('AvProj.mat')
catch
  load('AvProj_default.mat')
end
figure(1)
subplot(1,2,1)
plot(log(1:m_max-4), log(Proj_cst(1:m_max-4,40)),'r',...
  log(1:m_max-4), log(Proj_cst(1:m_max-4,35)),'b',...
  log(1:m_max-4), log(Proj_cst(1:m_max-4,30)),'g',...
  log(1:m_max-4), log(Proj_cst(1:m_max-4,25)),'k',...
  log(1:m_max-4), log(Proj_cst(1:m_max-4,20)),'y')
set(legend('n=40','n=35','n=30','n=25','n=20'),...
  'Location','northwest','Fontsize',16)
xlabel('Logarithm of m','Fontsize',16)
ylabel('Logarithm of \lambda','Fontsize',16)
title('Behavior of the projection constant \lambda wrto the subspace dimension m',...
  'Fontsize',14)
subplot(1,2,2)
plot(log(21:n_max), log(1/2-log(Proj_cst(16,21:n_max))/log(16)),'r',...
  log(21:n_max), log(1/2-log(Proj_cst(14,21:n_max))/log(14)),'b',...
  log(21:n_max), log(1/2-log(Proj_cst(12,21:n_max))/log(12)),'g',...
  log(21:n_max), log(1/2-log(Proj_cst(10,21:n_max))/log(10)),'k',...
  log(21:n_max), log(1/2-log(Proj_cst(8,21:n_max))/log(8)),'y')
set(legend('m=16','m=14','m=12','m=10','m=8'),...
  'Location', 'northeast','Fontsize',16)
xlabel('Logarithm of n','Fontsize',16)
ylabel('Logarithm of 1/2-ln(\lambda)/ln(m)','Fontsize',16)
title('Behavior of the projection constant \lambda wrto the superspace dimension n',...
  'Fontsize',14)


%% Empirical determination of projection constants in matrix spaces
%
% The code below exploits the function MinProjMat
% to compute the projection constants of certain subspaces of
% $\mathcal{M}_{\infty \to \infty}^{k \times k}$
% and then verifies the postulated behavior in $k$.

%% 
% We start with the subspace of symmetric matrices.

kMax = 10;
equalities = zeros(kMax,1);
tic;
for k = 2:kMax
  %define a basis for symmetric matrices as the columns of U
  U = zeros(k*k,k*(k+1)/2);
  for i = 1:k
    %express E_{i,i} as a column of U
    U(k*(i-1)+i,i) = 1;
  end
  for i = 1:k-1
    for j = i+1:k
      %express E_{i,j}+E_{j,i} as a column of U
      U(k*(i-1)+j,k+(k-1)*k/2-(k-i)*(k-i+1)/2+j-i) = 1;
      U(k*(j-1)+i,k+(k-1)*k/2-(k-i)*(k-i+1)/2+j-i) = 1;
    end
  end
  equalities(k) = ( abs(MinProjMatr(U,Inf) - (1+k)/2) < 1e-5 );
end
t = toc;
if sum(equalities) == kMax - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% The result is the same in $\mathcal{M}_{1 \to 1}^{k \times k}$

kMax = 10;
equalities = zeros(kMax,1);
tic;
for k = 2:kMax
  %define a basis for symmetric matrices as the columns of U
  U = zeros(k*k,k*(k+1)/2);
  for i = 1:k
    %express E_{i,i} as a column of U
    U(k*(i-1)+i,i) = 1;
  end
  for i = 1:k-1
    for j = i+1:k
      %express E_{i,j}+E_{j,i} as a column of U
      U(k*(i-1)+j,k+(k-1)*k/2-(k-i)*(k-i+1)/2+j-i) = 1;
      U(k*(j-1)+i,k+(k-1)*k/2-(k-i)*(k-i+1)/2+j-i) = 1;
    end
  end
  equalities(k) = ( abs(MinProjMatr(U,1) - (1+k)/2) < 1e-5 );
end
t = toc;
if sum(equalities) == kMax - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% Verification of the postulated behavior for the projection constant of the subspace
% of zero-trace matrices in $\mathcal{M}_{infty \to infty}^{k \times k}$

kMax = 10;
equalities = zeros(kMax-1,1);
tic;
for k = 2:kMax
  %define a basis for zero-trace matrices as the columns of U,
  %namely E_{i,j} for i~=j and E_{i,i}-E_{k,k} for i<k
  UU = zeros(k*k,k*k);
  for i = 1:k
    for j = 1:k
      UU(k*(i-1)+j,k*(i-1)+j) = 1;
      if i == j
        UU(k^2,k*(i-1)+j) = -1;
      end
    end
  end
  U = UU(:,1:end-1);
  equalities(k) = ( abs(MinProjMatr(U,Inf) - (2-1/k)) < 1e-5 );
end
t = toc;
if sum(equalities) == kMax-1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% The result is the same relative to $\mathcal{M}_{1 \to 1}^{k \times k}$

kMax = 10;
equalities = zeros(kMax-1,1);
tic;
for k = 2:kMax
  %define a basis for zero-trace matrices as the columns of U,
  %namely E_{i,j} for i~=j and E_{i,i}-E_{k,k} for i<k
  UU = zeros(k*k,k*k);
  for i = 1:k
    for j = 1:k
      UU(k*(i-1)+j,k*(i-1)+j) = 1;
      if i == j
        UU(k^2,k*(i-1)+j) = -1;
      end
    end
  end
  U = UU(:,1:end-1);
  equalities(k) = ( abs(MinProjMatr(U,1) - (2-1/k)) < 1e-5 );
end
t = toc;
if sum(equalities) == kMax-1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% Verification of the postulated behavior for the projection constant of the
% hyperplane orthogonal to shifted diagonal matrices in $\mathcal{M}_{infty \to infty}^{k \times k}$

kMax = 10;
equalities = zeros(kMax,kMax-1);
tic;
for k = 2:kMax
  for l = 1:k-1
    Dl = diag(ones(k-l,1),l);
    DlTranspose = Dl';
    DlVectorized = DlTranspose(:);
    [Q,R] = qr(DlVectorized);
    U = Q(:,2:end);
    equalities(k,l) = ( abs(MinProjMatr(U,Inf) - (2-1/(k-l))) < 1e-5 );
  end
end
t = toc;
if sum(sum(equalities)) == (kMax-1)*kMax/2
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% The result is the same relative to $\mathcal{M}_{1 \to 1}^{k \times k}$

kMax = 10;
equalities = zeros(kMax,kMax-1);
tic;
for k = 2:kMax
  for l = 1:k-1
    Dl = diag(ones(k-l,1),l);
    DlTranspose = Dl';
    DlVectorized = DlTranspose(:);
    [Q,R] = qr(DlVectorized);
    U = Q(:,2:end);
    equalities(k,l) = ( abs(MinProjMatr(U,1) - (2-1/(k-l))) < 1e-5 );
  end
end
t = toc;
if sum(sum(equalities)) == (kMax-1)*kMax/2
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%% 
% Verification of the postulated behavior for the projection constant of the
% hyperplane orthogonal to sums of rows of ones in $\mathcal{M}_{infty \to infty}^{k \times k}$

kMax = 10;
equalities = zeros(kMax,kMax);
tic;
for k = 2:kMax
  for l = 1:k
    sumRl = [ones(l,k); zeros(k-l,k)];
    sumRlTranspose = sumRl';
    sumRlVectorized = sumRlTranspose(:);
    [Q,R] = qr(sumRlVectorized);
    U = Q(:,2:end);
    equalities(k,l) = ( abs(MinProjMatr(U,Inf,'precision','low') - (2-2/l/k)) < 1e-5 );
  end
end
t = toc;
if sum(sum(equalities)) == kMax*(kMax+1)/2 - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% The result is slightly different for $\mathcal{M}_{1 \to 1}^{k \times k}$,
% as $max\{ 2,\ell \}$ come into play.

kMax = 10;
equalities = zeros(kMax,kMax);
tic;
for k = 2:kMax
  for l = 1:k
    sumRl = [ones(l,k); zeros(k-l,k)];
    sumRlTranspose = sumRl';
    sumRlVectorized = sumRlTranspose(:);
    [Q,R] = qr(sumRlVectorized);
    U = Q(:,2:end);
    equalities(k,l) = ( abs(MinProjMatr(U,1,'precision','low') - (2-2/max(2,l)/k)) < 1e-5 );
  end
end
t = toc;
if sum(sum(equalities)) == kMax*(kMax+1)/2 - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% Verification of the postulated behavior for the projection constant of the
% hyperplane orthogonal to sums of columns of ones in $\mathcal{M}_{infty
% \to infty}^{k \times k}$.
% Again, $max\{ 2,\ell \}$ come into play.

kMax = 10;
equalities = zeros(kMax,kMax);
tic;
for k = 2:kMax
  for l = 1:k
    sumCl = [ones(k,l) zeros(k,k-l)];
    sumClTranspose = sumCl';
    sumClVectorized = sumClTranspose(:);
    [Q,R] = qr(sumClVectorized);
    U = Q(:,2:end);
    equalities(k,l) = ( abs(MinProjMatr(U,Inf,'precision','low') - (2-2/max(2,l)/k)) < 1e-5 );
  end
end
t = toc;
if sum(sum(equalities)) == kMax*(kMax+1)/2 - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% The result is slightly different for $\mathcal{M}_{1 \to 1}^{k \times k}$,
% as max\{ 2,\ell \} does not come into play here.

kMax = 10;
equalities = zeros(kMax,kMax);
tic;
for k = 2:kMax
  for l = 1:k
    sumCl = [ones(k,l) zeros(k,k-l)];
    sumClTranspose = sumCl';
    sumClVectorized = sumClTranspose(:);
    [Q,R] = qr(sumClVectorized);
    U = Q(:,2:end);
    equalities(k,l) = ( abs(MinProjMatr(U,1,'precision','low') - (2-2/l/k)) < 1e-5 );
  end
end
t = toc;
if sum(sum(equalities)) == kMax*(kMax+1)/2 - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% Verification of the postulated behavior for the projection constant of the
% hyperplane orthogonal to the matrix with ones on the main diagonal and antidiagonal
% in $\mathcal{M}_{infty \to infty}^{k \times k}$

kMax = 10;
equalities = zeros(kMax,1);
tic;
for k = 2:kMax
  X = zeros(k,k);
  for i = 1:k
    X(i,i) = 1;
    X(i,k-i+1) = 1;
  end
  XTranspose = X';
  XVectorized = XTranspose(:);
  [Q,R] = qr(XVectorized);
  U = Q(:,2:end);
  equalities(k) = ( abs(MinProjMatr(U,Inf) - (2-1/k)) < 1e-5 );
end
t = toc;
if sum(sum(equalities)) == kMax - 1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%%
% The result is the same relative to $\mathcal{M}_{1 \to 1}^{k \times k}$

kMax = 10;
equalities = zeros(kMax,1);
tic;
for k = 2:kMax
  X = zeros(k,k);
  for i = 1:k
    X(i,i) = 1;
    X(i,k-i+1) = 1;
  end
  XTranspose = X';
  XVectorized = XTranspose(:);
  [Q,R] = qr(XVectorized);
  U = Q(:,2:end);
  equalities(k) = ( abs(MinProjMatr(U,1) - (2-1/k)) < 1e-5 );
end
t = toc;
if sum(sum(equalities)) == kMax-1
  disp(strcat('Values of projection constants verified up to k=', ...
    num2str(kMax), 32, '(', num2str(round(t/60)), 'min)'))
end

%% Experiments in polynomial spaces

%%
% Determination of the projection constant for the 2-codimensional subspace
% of polynomials vanishing at the endpoints

x = chebfun('x');
nMin = 3; nMax = 11;
nu = 200;
lowerBounds = zeros(1,nMax);
upperBounds = zeros(1,nMax);
tic;
for n = nMin:nMax
  m = n-2;
  basis = cell(1,m);
  for j = 1:m
    basis{j} = (1-x.^2).*x.^(j-1);
  end
  bounds = MinProjPoly(basis,n,nu,'SDP');
  lowerBounds(n) = bounds(1);
  upperBounds(n) = bounds(2);
end
t = toc;

%%
% Determination of the projection constant of the space P^m as a subspace of P^n

x = chebfun('x');
nMin = 3; nMax = 11;
nu = 200;
lowerBounds = zeros(nMax,nMax-1);
upperBounds = zeros(nMax,nMax-1);
tic;
for n = nMin:nMax
  for m = 2:n-1
    basis = cell(1,m);
    for j = 1:m
      basis{j} = x.^(j-1);
    end
    bounds = MinProjPoly(basis,n,nu,'SDP');
    lowerBounds(n,m) = bounds(1);
    upperBounds(n,m) = bounds(2);
  end
end
t = toc;

%%
% Verification that the minimal projection from P^(m+1) to P^(m) is interpolating

x = chebfun('x');
m = 5;
n = m+1;
nu = 300;
basis = cell(1,m);
for j = 1:m
  basis{j} = x.^(j-1);
end
tic;
[bounds1,P1] = MinProjPoly(basis,n,nu,'LP');
[bounds2,P2] = MinProjPoly(basis,n,nu,'SDP');
t = toc;
P1Tm = chebfun(P1*(-1).^((n+1:2*n)'));
fullP2 = full(P2);
P2Tm = chebfun(fullP2(:,n),'coeffs');
figure(2);
subplot(1,2,1); plot(P1Tm-chebpoly(m),'b');
set(legend('P^{LP}(T_m)-T_m'),'Location','southeast','Fontsize',16)
title(strcat(num2str(bounds1(1)),'<||P^{LP}||<',num2str(bounds1(2))),'Fontsize',16)
subplot(1,2,2); plot(P2Tm-chebpoly(m),'r');
set(legend('P^{SDP}(T_m)-T_m'),'Location','southwest','Fontsize',16)
title(strcat(num2str(bounds2(1)),'<||P^{SDP}||<',num2str(bounds2(2))),'Fontsize',16)

%% Computation of the minimal extension of a specific map onto polynomial spaces

%%
% The map defined by the Legendre extrema when l=3
l = 3;
nMin = 4; nMax = 14;
extCst = zeros(1,nMax);
tic;
for n = nMin:nMax
  pts = [-1; jacpts(n-2,1,1); 1];
  extCst(n) = ExMinExt(l,pts);
end
t = toc;

%%
% The map defined by the Chebyshev extrema when l=4
l = 4;
nMin = 4; nMax = 14;
extCst = zeros(1,nMax);
tic;
for n = nMin:nMax
  pts = chebpts(n);
  extCst(n) = ExMinExt(l,pts);
end
t = toc;

%% References
%
% 1. S. Foucart, "Computation of minimal projections and extensions",
% In preparation.
%
% 2. CVX Research, Inc., "CVX: MATLAB software for disciplined
% convex programming, version 2.1", 2014, http://cvxr.com/cvx.
%
% 3. L. N. Trefethen et al., "Chebfun Version 5, The Chebfun Development Team", 2014,
% http://www.chebfun.org.
%
% 4. J. Blatter and E. W. Cheney, "Minimal projections on hyperplanes in
% sequence spaces", Annali di Matematica Pura ed Applicata 101.1 (1974): 215-227.
%
% 5. G. Lewicki, "Minimal Projections onto Two Dimensional Subspaces of $\ell_\infty^{(4)}$",
% Journal of Approximation Theory, 88.1 (1997): 92-108.
