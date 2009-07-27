%
% Copyright Alan C. Evans
% Professor of Neurology
% McGill University
%
function length = laplace_pathlength(ly, lx)

S = max(size(ly));
length = 0;
for i=1:S-1
  % hello pythagoras
  length = length + sqrt( ( abs(ly(i+1) - ly(i)))^2 + ( abs(lx(i+1) - ...
                                                 lx(i)))^2);
end
