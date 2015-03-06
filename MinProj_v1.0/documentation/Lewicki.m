%%
% Lewicki.m
% Computes the projection constant of 2-dimensional subspace of
% $\ell_\infty^4$ according to the formula given in
% MINIMAL PROJECTIONS ONTO TWO DIMENSIONAL SUBSPACES OF \ELL_\INFTY^4
% by Lewicki

% Find solutions to
% minimize ||P||_{infty -> infty}
% over all projections P from R^4 onto a given 2-dimensional subspace 
%
% Usage: [projCst, isApplicable] = Lewicki(f,g,verb)
%
% f,g: vectors in R^4 orthogonal to the subspace under consideration
% verb: 0 or 1 if one doesn't/does want a message displaying the failing condition
% 
% projCst: the projection constant of the subspace relative to the infinity-norm
% isApplicable: 1 if the conditions of applicability of the formula hold,
%               0 if they do not hold
%
% Written by Simon Foucart in June 2014
% Last updated in December 2014
% Send comments to simon.foucart@centraliens.net

function [projCst, isApplicable] = Lewicki(f,g,verb)

if nargin<3
    verb = 1;  % by default, a message is displayed
end

% normalize f and g
f=f/norm(f,1); g=g/norm(g,1);

% verify the conditions of applicability of the formula hold
if f(1)<=0 || f(2)~=0 || f(3)<=0 || f(4)<=0
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (2.1) fails for f')
    end
    return
end

if g(1)~=0 || g(2)<=0 || g(3)<=0 || g(4)<=0
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (2.1) fails for g')
    end
    return
end

if f(1)>=1/2 || f(3)>=1/2 || f(4)>=1/2
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (2.2) fails for f')
    end
    return
end

if g(2)>=1/2 || g(3)>=1/2 || g(4)>=1/2
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (2.2) fails for g')
    end
    return
end
    
if abs(f(1)*g(3)-f(3)*g(2)) >= abs(f(3)*g(4)-f(4)*g(3))
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (2.3) fails')
    end
    return
end

if abs(f(1)*g(4)-g(2)*f(4)) >= abs(f(3)*g(4)-f(4)*g(3))
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (2.4) fails')
    end
    return
end

if f(3)<=g(3) || f(3)>g(2)+g(3)*(1-2*f(1)) || g(4)>f(1)+f(4)*(1-2*g(2))
    projCst = NaN;
    isApplicable = 0;
    if verb == 1
        disp('Error: condition (3.1) fails')
    end
    return
end

% when the conditions hold, write the formula for the projection constant

a = 1 + 1/...
    ( g(2)/(1-2*g(2)) + g(3)*(1-2*f(1))/(1-2*g(2))/(1-2*f(3)) + g(4)/(1-2*g(4)) );
b = 1 + 1/...
    ( f(1)/(1-2*f(1)) + f(4)*(1-2*g(2))/(1-2*f(1))/(1-2*g(4)) + f(3)/(1-2*f(3)) );

projCst = max(a,b);
isApplicable = 1;
 
end