function dTheta_dt = pendulum_linearized_1(t, Theta, A_big)
    %computes the derivative of the state vector for a linearized N-link pendulum.
    dTheta_dt = A_big * Theta;
end