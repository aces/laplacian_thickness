%
% Copyright Alan C. Evans
% Professor of Neurology
% McGill University
%
function thickness_lattice = laplace_thickness(lattice, dx, dy)
% LAPLACE_THICKNESS - create a laplacian thickness map in 2 dimensions
% 
% [THICKNESS_LATTICE] = LAPLACE_THICKNESS(LATTICE, DX, DY)
% Evaluates the length of the streamline along the gradients DX and
% DY on the LATTICE. The gradients are generated using
% LAPLACE_GRADIENT and the lattice is created by SOLVE_LAPLACE2

% select all points to solve over
[Fy,Fx] = find(lattice<10000 & lattice>0);

% initialise output
thickness_lattice = zeros(size(lattice));
size(thickness_lattice)
for i=1:length(Fx)
  %integrate towards 0
  [lx,ly] = laplace_euler(Fy(i),Fx(i),10000,dy,dx,1,lattice);
  length = laplace_pathlength(ly,lx);

  %integrate towards 10000
  [lx,ly] = laplace_euler(Fy(i),Fx(i),0,dy,dx,-1,lattice);
  length = length + laplace_pathlength(ly,lx)
  i
  thickness_lattice(Fy(i),Fx(i)) = length;
end
