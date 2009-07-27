%
% Copyright Alan C. Evans
% Professor of Neurology
% McGill University
%
function [dx, dy] = laplace_gradient(lattice)
% LAPLACE_GRADIENT - gradient of laplace grid
%
% [DX, DY] = LAPLACE_GRADIENT(LATTICE)
% Evaluates the gradient of the LATTICE, returning it in the X and
% Y dimension using the two point difference formula. The LATTICE
% is generated using SOLVE_LAPLACE2

%use two point difference formula

S = size(lattice);

dx = ones(S);
dy = ones(S);



for i=2:S(1)-1
  for j=2:S(2)-1
    dy(i,j) = (lattice(i+1,j) - lattice(i-1,j)) / 2;
    dx(i,j) = (lattice(i,j+1) - lattice(i,j-1)) / 2;
  end
end

% normalise 
dx = dx ./ sqrt(dx .^ 2 + dy .^ 2);
dy = dy ./ sqrt(dy .^ 2 + dx .^ 2);
