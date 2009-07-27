%
% Copyright Alan C. Evans
% Professor of Neurology
% McGill University
%
function grid = laplace_prepare_grid(S, S_prime)
% LAPLACE_PREPARE_GRID - initialise grid for solving
% GRID = LAPLACE_PREPARE_GRID(S, S_prime)
%   GRID = the initialised grid, where values to be solved are
%          assigned a value of 5000, the outside boundary a value of
%          10000, and the inside boundary a value of 0.
%   S = vector of integers defining one of the boundaries
%   S_prime = vector of integers defining the other boundary.
%
% See also SOLVE_LAPLACE2

size1 = size(S);
size2 = size(S_prime);

if size1 ~= size2
  error('Surfaces not the same size');
end

v_size = (max(S) - min(S_prime)) + 15;
h_size = size1(2);

S = (S - min(S_prime)) + 8;
S_prime = (S_prime - min(S_prime)) + 8;

grid = ones(v_size,h_size);
grid = grid + 4999;
for i=1:h_size
  grid(S(i)+1:v_size,i) = 10000;
  grid(1:S_prime(i)-1,i) = 0;
end


