function dTheta_dt = pendulum_ct(t, Theta, n, g, m, l, C, omega_0, l_max, l_min, conc, ct_span)

    theta = Theta(1:n);
    omega = Theta(n+1:2*n);
    F_ext = C * cos(omega_0 * t);
    l_dot = zeros(n);
    l_2dot = zeros(n);
    
    % change the last length as a function of t
    xL = -0.02:.002:0.02;
    xR = 0.98:.002:1.02;
    order = 2;
    tol = 3e-3;
    concavity = conc;
    slope = 1;
    xL1 = -0.05:.05:0.05;

    x = [xL, 0.5-concavity+xL1, xR];
    y = [zeros(size(xL)), (0.5+concavity)+xL1*slope, ones(size(xR))];

    w = ones(size(x));
    w([1:length(xL),length(x)-length(xR)+1:length(x)]) = 10;
    sp = spaps(x,y, tol, w, order);

    % 1st & 2nd order diff of the last length
    if ct_span(1)<t && t < ct_span(2)
        l(n) = (l_max-l_min) * fnval(sp, (t-ct_span(1))/(ct_span(2)-ct_span(1))) + l_min;
        f1 = fnder(sp);
        f2 = fnder(sp,2);
        l_dot(n) = (l_max-l_min) * fnval(f1, (t-ct_span(1))/(ct_span(2)-ct_span(1)))/(ct_span(2)-ct_span(1));
        l_2dot(n) = (l_max-l_min) * fnval(f2, (t-ct_span(1))/(ct_span(2)-ct_span(1)))/(ct_span(2)-ct_span(1))^2;
    end
    if t<ct_span(1)
        l(n) = l_min;
        l_dot(n) = 0;
        l_2dot(n) = 0;
    end
    if t>ct_span(2)
        l(n) = l_max;
        l_dot(n) = 0;
        l_2dot(n) = 0;
    end
    
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
            first_term_of_B_q = first_term_of_B_q + ...
                (mass_sum * ((l(q) * l(k) * omega(k)^2 +l(q) * l_2dot(k))* sin(theta(q) - theta(k))...
                + 2 * omega(k) * l_dot(k) * l(q) * cos(theta(q)-theta(k))));
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