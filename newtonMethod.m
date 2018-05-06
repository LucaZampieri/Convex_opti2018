function [xOpt, i] = newtonMethod(f, grad, hess, x, tol, maxIter)

xOpt = x; % initial solution

% Parameters for backtracking
alpha = 0.01; % alpha in (0, 0.5)
beta = 0.5; % beta in (0, 1)

for i = 1:maxIter
    
    dx = -hess(xOpt)\grad(xOpt); % Compute increment
    lambda = -grad(xOpt)'*dx;
    
    if lambda/2 <= tol
        
        fprintf('Newton method convergence reached in %d iteration\n', i);
        return
        
    else
        
        t = backtrackingLineSearch(f, grad, xOpt, dx, alpha, beta);
        xOpt = xOpt + t*dx; % Update solution
        
    end
    
end

fprintf('Newton method convergence not reached after %d iteration\n', i);

end

