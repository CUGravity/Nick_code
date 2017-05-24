function testUKF
tic
increment = 0.01; %1e14*eps;
[x_kp, p, w] = initialize(increment); % initial x0 et alia
P_km = p.P;
x_km = x_kp + [diag(P_km(1:3,1:3));0; diag(P_km(4:6,4:6))]*0.5;
orbits = 1; % number of desired orbits to simulate
time = 0;   tend = 4; %1200; %3600; %3600*orbits; 
counter = 1;



    while time < tend
        time_prev = time
        time = time + increment;
        
        %% get measurements 
        x_measured =  getGyroData(x_kp, time);
        
        %% compute kalman filter update
        [x_koutput, P_k, K_k] = UKF(x_km, P_km, w, x_measured, p.H, p.R, p.Q, increment);
        P_km = P_k;
        x_km = x_koutput;
        x_koutput = normalize(x_koutput);
        x_kalman(counter,:) = x_koutput';
        covariance(counter,:) = diag(P_km);
        gain(counter,:) = diag(K_k);
        
        x_k = updatedSYS(time, time_prev, x_kp, w);
        x_k = normalize(x_k);
        x_system(counter,:) = x_k';
        t(counter) = time;
        
        %% update for next iteration
        x_kp = x_k;
        counter = counter+1;
    end
    
    clf; 
    figure(1)
    subplot(2,1,1)
    plot(t,x_system(:,1),t,x_system(:,2),t,x_system(:,3),t,x_system(:,4));
    title('quat')
    legend('q_1', 'q_2', 'q_3', 'q_4')
    
    subplot(2,1,2)
    plot(t,x_system(:,5),t,x_system(:,6),t,x_system(:,7));
    title('bias')
    legend('b_1', 'b_2', 'b_3')
    
    figure(2)
    subplot(2,1,1)
    plot(t,x_kalman(:,1),t,x_kalman(:,2),t,x_kalman(:,3),t,x_kalman(:,4));
    title('quat_k')
    legend('q_1k', 'q_2k', 'q_3k', 'q_4k')
    
    subplot(2,1,2)
    plot(t,x_kalman(:,5),t,x_kalman(:,6),t,x_kalman(:,7));
    title('bias_k')
    legend('b_1k','b_2k','b_3k')
      
    finaltime = toc/60
end



function x_normed = normalize(x)
%% Normalize quaternion
q = x(1:4);
x(1:4) = q/norm(q);
x_normed = x;
end

function x_measured = getGyroData(x_kp, t_k)
%% Spoof sensors to get measurements
% state array is: x, xdot, w, q (3+3+3+4)x1 = 13x1
%% FOR NOW:
gyroNoise = 0.875;
quatNoise = 0.01;
% x_measured = x_kp + [randn(3,1)*gyroNoise; randn(4,1)*quatNoise];
x_measured = x_kp(1:4);
end

function x_k = updatedSYS(t_k, t_km, x_km, w)
%% Propagate dynamics to state array
    eta_1 = 1.e-10;
    eta_2 = 1.e-7;

    dt = t_k-t_km;
    dx = SYSdynamics(x_km, w, dt); 
    x_k = x_km + dx*dt;
    x_k(5:end) = x_k(5:end) - randn(3,1)*eta_2; % normal variable w/ mu = 0 and stdev = eta_v
    w = w - x_k(5:end) - randn(3,1)*eta_1;
    
    %% Re-normalize q
    x_k = normalize(x_k);
%     x_k = dx;
end

function dx =  SYSdynamics(x_km, w, increment)
n = length(x_km);
dx = zeros(size(x_km));
    q1 = x_km(1);
    q2 = x_km(2);
    q3 = x_km(3);
    q4 = x_km(4);

%% Compute update to state array
    q = [q1 q2 q3 q4]';
    qdot = 0.5*[crs(w) w; -w' 0]*q;
end

function [x0, p, w0] = initialize(dt)
    a = [1; -2; 1]/norm([1; -2; 1]);
    phi = 0;
    q0 = [a*sin(phi/2); cos(phi/2)];
    w0 = a*0.4;
    b0 = [0;0;0];
    x0 = [q0; b0];
    
    sigma_u = 1.e-10; % white noise stdev (eta_1)
    sigma_v = 1.e-7; % bias stdev (eta_2)
    meas_sigma = 9e-4; % error with measurement sensors
    
    p.H = eye(6);  % sensor to state transformation matrix 
    p.R = eye(4)*meas_sigma^2;  % Sensor error covariance matrix
    p.P = 2*diag([0.009*ones(3,1); (9.7e-7)*ones(3,1)]);  % initial Process covariance matrix estimate
    p.Q = 1/2*dt*diag([ones(3,1)*(sigma_v^2-1/6*dt^2*sigma_u^2); sigma_u^2*ones(3,1)]);  % process error covariance matrix
    % 1e14*
end