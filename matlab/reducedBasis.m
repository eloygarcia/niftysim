function [Phi,lam] = reducedBasis(U,M)
% REDUCEDBASIS - Compute the M-dimensional reduced basis for the
% displacement matrix U. U (3N-by-S) is an ensemble of S displacement
% vectors, each of length 3N. The displacements would normally be sampled
% from the complete displacement history of a model with N nodes.
%
% Usage:
% [Phi,lam] = reducedBasis(U,M)
%
% Inputs:
% U: matrix of nodal displacements at various time points (3N-by-S)
% M: dimensionality of the reduced basis (scalar)
%
% Returns:
% Phi: matrix of eigenvectors constituting the reduced basis (3N-by-M)
% lam: vector of eigenvalues corresponding to the eigenvectors in Phi
%      (3N-by-1)
%
% Copyright (c) 2010, University of Queensland. All rights reserved.
% MedTeQ Centre
% See the LICENSE.txt file in the root folder
% ztaylor@itee.uq.edu.au
% http://www.itee.uq.edu.au/~medteq

method = 1;

% Compute row means
disp('Compute row means')
ubar = mean(U,2);

% Centre the U matrix
disp('Centre the U matrix')
[N,S] = size(U);
% Ubar = zeros(N,S);
for i = 1:S
   U(:,i) = U(:,i) - ubar;
end

if N < S % Direct method
   % Compute covariance matrix (left singular)
   disp('Compute covariance matrix')
   Cd = U*U'/M;

   % Solve the eigenproblem
   disp('Solve the eigenproblem')
   [Phi,D] = eigs(Cd,M);
   lam = diag(D);
else % Method of snapshots - more common case
   % Compute covariance matrix (right singular)
   disp('Compute covariance matrix')
   Cd = U'*U;

   % Solve the eigenproblem
   disp('Solve the eigenproblem')
   [Phi,D] = eigs(Cd,M);
   lam = diag(D);
   % Scale the eigenvectors
   for i = 1:M
      Phi(:,i) = Phi(:,i)/sqrt(lam(i));
   end
   % Project onto displacements
   Phi = U*Phi;
   % Scale eigenvalues
   lam = lam/M;
end