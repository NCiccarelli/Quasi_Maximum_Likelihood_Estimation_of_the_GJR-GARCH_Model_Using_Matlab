function [c,ceq] = mycon_threshold_garch(x)


c(1) = ( (x(1)) + (x(2))  ) .*0.5 .*  (x(4).^2) +   x(3)   - 0.999999999999999 ; % Nonlinear inequalities at x.

 ceq = [] ;    % Compute nonlinear equalities at x.
end
 

 
 
 
 
 
 
 
 
