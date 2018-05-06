function [xOpt, newtonIt] = barrierMethod(f0, grad_f0, hess_f0, phi, grad_phi, hess_phi, x, t, mu, m, tol, maxIter, tolNewton, maxIterNewton)


for i = 1:maxIter
    
    f = @(x) t * f0(x) + phi(x);
    grad = @(x) t * grad_f0(x) + grad_phi(x);
    hess = @(x) t * hess_f0(x) + hess_phi(x);
    
    [xOpt, newtonIt(i)] = newtonMethod(f, grad, hess, x, tolNewton, maxIterNewton);
    x = xOpt;
    
    if m/t < tol
        
        fprintf('Barrier method convergence reached in %d iteration\n', i);
        return
        
    else
        
        t = mu * t;
        
    end
    
    
end

fprintf('Barrier method convergence not reached after %d iteration\n', i);

end


