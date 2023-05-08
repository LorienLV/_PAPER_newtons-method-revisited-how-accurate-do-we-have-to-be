function [x, rel]=newton_sqrt(s,maxit,seed,epsilon)

% NEWTON_SQRT  Demonstrates the resilience of Newtons method
%
% CALL SEQUENCE: [x, rel]=newton_sqrt(s,maxit,seed,delta)
%
% INPUT:
%   s      a vector containing real numbers in [1,4]
%   maxit  do maxit Newton steps
%   seed   a seed for the random generator
%   delta  the relative errors are at least 0.5*delta and at most delta
%
% OUTPUT:
%   x     a column vector containing square roots of s
%   rel   a matrix such that 
%           rel(i,j+1) is the relative error for sqrt(s(i))
%           after j Newton steps
%
% MINIMAL WORKING EXAMPLE: figure1.m

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2021-07-12  Initial programming and testing
%   2023-05-08  Updated the inline documentation

% Seed the random number generator
rng(seed);

% Isolate number of elements
n=numel(s);

% Reshape as column vector
s=reshape(s,n,1);

% Allocate space for all errors
rel=zeros(n,maxit+1);

% Compute a good initial guess for sqrt(s), s in [1,4]
% This approximation is the best uniform linear approximation
x=1/3*s+17/24;

% Target values
z=sqrt(s);

% Initial relative error
rel(:,1)=(z-x)./z;

% Main loop
for j=1:maxit
    % Compute the correct correction
    dx=-(x.^2-s)./(2*x);
        
    % Make random numbers in the interval [1/2;1]
    aux=0.5*(1+rand(n,1)); 
    
    % Flip some signs
    aux=aux.*sign(rand(n,1)-rand(n,1)); 
   
    % Deliberately contaminate the correction
    dx=dx.*(1+epsilon*aux);
    
    % Compute the next approximation, 
    % but use the contaminated correction
    x=x+dx;
    
    % Compute the relative error
    rel(:,j+1)=(z-x)./z;
end
    