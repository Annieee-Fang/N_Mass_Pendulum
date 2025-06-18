function dTheta_dt = pendulum_1(t, Theta, n, g, m, l, C, omega_0)

    theta = Theta(1:n);
    omega = Theta(n+1:2*n);
    F_ext = C * cos(omega_0 * t);
    
    A = zeros(n, n);
    B = zeros(n, 1);
    
    for q = 1:n

        for k = 1:n
            mass_sum = sum(m(max(q, k):n)); 
            A(q, k) = mass_sum * l(q) * l(k) * cos(theta(q) - theta(k));
        end
        
        first_term_of_B_q = 0;
        for k = 1:n
            mass_sum = sum(m(max(q, k):n));
            first_term_of_B_q = first_term_of_B_q + (mass_sum * l(q) * l(k) * omega(k)^2 * sin(theta(q) - theta(k)));
        end
        if q == 1
                first_term_of_B_q = first_term_of_B_q + F_ext;
        end
        
        second_term_of_B_q = g * l(q) * sin(theta(q)) * sum(m(q:n));
        
        B(q) = first_term_of_B_q + second_term_of_B_q;
    end
    
    omega_dot = A \ -B;
    
    dTheta_dt = [omega; omega_dot];
end