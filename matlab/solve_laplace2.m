function lattice = solve_laplace2(grid, convergence_criteria)
% SOLVE_LAPLACE2 - solves laplaces equation in two dimensions
% LATTICE = SOLVE_LAPLACE2(S, S_PRIME, CONVERGENCE_CRITERIA)
% SOLVE_LAPLACE2 takes two two-dimensional vectors as an argument,
% each representing the x,y coordinates of a surface. It then
% solves laplace's equation between these surfaces, returning the
% lattice. The last argument is the convergence 


%size(grid)
%for i=1:size(S)
%  points(S_prime(i,1):S(i,1),i) = 1;
%  % and intialise the grid to a reasonable guess
%  f = find(grid(:,i)==1);
%  num_points = size(f);
%  scale = S_bound / num_points(1);
%  for j=1:num_points
%    grid(f(j),i) = j*scale;
%  end
  
%end
%points


%for i=1:size(S)
%  grid(S(i,1),S(i,2)) = S_bound;
%end
%for i=1:size(S_prime)
%  grid(S_prime(i,1),S_prime(i,2)) = S_prime_bound;
%end



S_bound = 10000;
S_prime_bound = 0;
grid2 = grid;
points = grid;
con_old = 10;
con = 10;
% now approximate laplaces equation
test = Inf;
grid_size = size(grid);
l = 0

%run this loop while the convergence criteria is not met
while test > convergence_criteria | test == Inf | test == NaN
  l = l+1  
  for i=2:grid_size(1)-1
    for j=2:grid_size(2)-1
      if points(i,j) == 5000 % only check the points we want
        if grid(i,j) == S_bound | grid(i,j) == S_prime_bound
          % keep points that have boundary values
          grid2(i,j) = grid(i,j);
        else
          % and, finally, solve the equation at that location
          grid2(i,j) = (grid(i+1,j) + grid(i,j+1) + grid(i-1,j) + ...
                        grid(i,j-1)) / 4;
        end
      end
    end
  end
  % check for convergence
  test = sum(sum(abs(grid2 - grid))) 

  grid = grid2;
end

lattice = grid;

