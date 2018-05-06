function t = backtrackingLineSearch(f, grad, x, dx, alpha, beta)

t = 1;

% Check that if new solution is a real number (e.g we want positive
% argument for the log)
while isreal(f(x+t*dx)) == false
    
    t = beta*t;
    
end


while f(x+t*dx) > f(x) + alpha * t * grad(x)' * dx;
    
    t = beta*t;
    
end





