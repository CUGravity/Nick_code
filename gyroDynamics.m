function [chi_kp,  q_kp] = gyroDynamics(chi, w, q_kp_0, q_km_0, increment)
%% Gyro DYNAMICS OF THE SYSTEM for CU Gravity EKF
% Equations taken in part from "Kalman Filtering for Spacecraft Systems 
% Alignment Calibration" by Mark E. Pittelkau from the Journal of
% Guidance, Control, and Dynamics (Vol. 24, No. 6, November-December
% 2001), and "Kalman Filtering fro Spacecraft Attitude Estimation" by 
% Lefferts, Markley, and Shuster from the Journal of Guidance, Control, 
% and Dynamics (Vol. 5, No. 5, September-October 1982).
%
%
% Contact Niccolo Porcari, npp26@cornell.edu, MEng 2016 Cornell University
    a = 1;
    f = 2*(a+1);
    %% Compute error quaternion
    delta_q4 = -a*norm(chi(1:3))^2+f*sqrt(f^2+(1-a^2)*norm(chi(1:3))^2)...
        /(f^2+norm(chi(1:3))^2);
    delta_q = 1/f*(a+delta_q4)*chi(1:3);
    detlaq_k = [delta_q; delta_q4];
    
    %% Compure error quaternion change
    q_k = [[delta_q4*eye(3) - crs(delta_q); -delta_q'] detlaq_k]*q_km_0;
    q_k = q_k/norm(q_k);
    
    %% Update error quaternion result 
    w = w + chi(4:end); % update omega with bias
    psi = sin(0.5*increment*norm(w))/norm(w)*w;
    q_kp = [cos(0.5*increment*norm(w))*eye(3)-crs(psi) psi;
            -psi' cos(0.5*increment*norm(w))]*q_k;
    q_kp = q_kp/norm(q_kp);
    deltaq_kp = [[q_kp_0(4)*eye(3) - crs(q_kp_0(1:3)); -q_kp_0(1:3)'] q_kp_0]...
        *[-q_kp(1:3); q_kp(4)];
    
%     deltaq_kp = [[q_kp(4)*eye(3) - crs(q_kp(1:3)); -q_kp(1:3)'] q_kp]...
%         *[-q_kp_0(1:3); q_kp_0(4)];
    
    chi_kp =[f * deltaq_kp(1:3)/(a+deltaq_kp(4)); chi(4:end)];
    
end