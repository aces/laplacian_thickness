function [lx,ly] = laplace_euler(y0,x0,stop_criteria,dy,dx,h,lattice)
% LAPLACE_EULER - create a streamline to one surface
% [LX,LY] = LAPLACE_EULER(Y0,X0,STOP_CRITERIA,DY,DX,H,LATTICE)
%
% This function returns two vectors, LX AND LY, defining the path
% taken to reach the the value STOP_CRITERIA in the potential field
% LATTICE from the starting points Y0 and X0 using the step size H.
%
% See also LAPLACE_THICKNESS, SOLVE_LAPLACE2, LAPLACE_GRADIENT

lx(1) = x0;
ly(1) = y0;
i = 1;
test_condition = stop_criteria +1;
s = size(lattice);
while test_condition ~= stop_criteria
%  ly
%  lx

  lx(i+1) = lx(i) + interp2(dx,lx(i),ly(i)) * h;
  ly(i+1) = ly(i) + interp2(dy,lx(i),ly(i)) * h;

  if lx(i+1) < 1
    lx(i+1) = 1;
  end
  if ly(i+1) < 1
    ly(i+1) = 1;
  end
  
 
  i = i+1;

  % abort the search if we have a NaN, else check to see wether
  % we've reached the stop criteria
  if (isnan(ly(i)) | isnan(lx(i)) | round(ly(i)) > s(1) | round(lx(i)) > s(2))
    [ly(i) lx(i)]
    test_condition = stop_criteria;
  else
    test_condition = lattice(round(ly(i)),round(lx(i)));
  end
end
