function [x_k, P_k, K_k] = UKF(x_km, P_km, w, y_k, H_k, R_k, Q_km, dt)
% Unscented Kalman Filter Code
% Contact: Niccolo Porcari 2017
    %% Needed Matrices
%     H_k = p.H;
%     R_k = p.R;  % Sensor error covariance matrix
%     Q_km= p.Q;  % process error covariance matrix
    %% Set parameters
    n = 6;
    lambda = 1;    
    a = 1;
    f = 2*(a+1);
    
    %% Compute Sigma Points
    sigma_k = sqrt(n+lambda)*chol((P_km+Q_km));
    chi_0 = [0;0;0; x_km(5:end)];
    xi_km = [chi_0, chi_0+sigma_k, chi_0-sigma_k];
    xi_k(:,1) = [zeros(3,1); chi_0(4:end)];
    
    w = w - xi_k(4:end); % update omega with bias
    psi = sin(0.5*dt*norm(w))/norm(w)*w;
%     x_km(1:4)
    q_kp_0 = [cos(0.5*dt*norm(w))*eye(3)-crs(psi) psi;
            -psi' cos(0.5*dt*norm(w))]*x_km(1:4);
    
    %% Propagate Sigma Points    
    for i = 2:2*n+1
        [xi_k(:,i), q_kp(:,i)]  = gyroDynamics(xi_km(:,i), w, q_kp_0, x_km(1:4), dt);
    end
    
    x_hat_kp = (lambda*xi_k(:,1)+1/2*sum(xi_k(:,2:end),2))./(n+lambda);
%     q_kp
%     xi_k
    %% Compute Necessary matrices for Kalman Gain
    SUM_term = zeros(6,6);
    for i = 2:2*n+1
       SUM_term = SUM_term + (xi_k(:,i)-x_hat_kp)*(xi_k(:,i)-x_hat_kp)';
    end
%     Q_km
    P_kp = (lambda*(xi_k(:,1)-x_hat_kp)*(xi_k(:,1)-x_hat_kp)' + ...
        1/2*SUM_term)./(n+lambda) + Q_km;
    
%     q_kp
    y_hat_kp = (lambda*q_kp_0+1/2*sum(q_kp,2))./(n+lambda);     % mean observation
    
    SUM_term = zeros(4,4);
    for i = 2:2*n+1
       SUM_term = SUM_term + (q_kp(:,i)-y_hat_kp)*(q_kp(:,i)-y_hat_kp)';
    end

    P_yy = (lambda*(q_kp_0-y_hat_kp)*(q_kp_0-y_hat_kp)' + ...
        1/2*sum(SUM_term,2))./(n+lambda);
    
    P_vv = P_yy + R_k; % here Rk is a 4x4
     
    SUM_term = zeros(6,4);
    for i = 2:2*n+1
       SUM_term = SUM_term + (xi_k(:,i)-x_hat_kp)*(q_kp(:,i)-y_hat_kp)';
    end
    
    P_xy = (lambda*(xi_k(:,1)-x_hat_kp)*(q_kp_0-y_hat_kp)' + ...
        1/2*sum(SUM_term,2))./(n+lambda);
    
    %% Gain Computation
    K_k = P_xy  * inv(P_vv);
    diff_k = y_k - y_hat_kp;     % measurement residual (i.e. error of prediction)

    %% Update estimate and estimate covariance
    q_x_k = x_hat_kp + K_k*diff_k;
%     P_kp
    P_k = P_kp - K_k*P_vv*transpose(K_k);
    
    %% Compute actual propagating state
    delta_q4_kp = -a*norm(q_x_k(1:3))^2 + f*sqrt(f^2+(1-a^2)*norm(q_x_k(1:3))^2)...
        /(f^2+norm(q_x_k(1:3))^2);
    delta_rho_kp = 1/f*(a+delta_q4_kp) * q_x_k(1:3);
    deltaq_kp = [delta_rho_kp; delta_q4_kp];
    x_k = [[deltaq_kp(4)*eye(3) - crs(deltaq_kp(1:3)); -deltaq_kp(1:3)'] deltaq_kp]*q_kp_0;
    x_k = x_k/norm(x_k);
    x_k = [x_k; q_x_k(4:end)];

end