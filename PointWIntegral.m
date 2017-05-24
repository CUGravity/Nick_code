function PointWIntegral
    clf; tic;
    %% desired is 4 RPM which translates to
    rate = 0.4; %4*2*pi/60; %
    axis_0 = [-0.2;-0.1;0.3]; %[1;0;0];%
    initial_spin = 0.05;
    
    a = axis_0/norm(axis_0);
    phi = 0.1;
    q0 = [a*sin(phi/2); cos(phi/2)];
    w0 = a*initial_spin;
    
    wc_N = [1;0;0]*rate;
    q_match = [[1;0;0]*sin(phi/2); cos(phi/2)];
    xi0 = q0-q_match;
    xi = xi0;
    x_kp = [w0;
            q0]; % 7x1 vector x2 (state+integral error)
   
    Inertia = diag([1500; 1200; 1200]);
    Inertia = Inertia; % + [0 10 -10; 0 10 0; -10 0 0];
    
    orbits = 1; % number of desired orbits to simulate
    time = 0;   tend = 400; %1200; %3600; %3600*orbits; 
    counter = 1;
    increment = 0.01; %1e14*eps;

    while time < tend
        time_prev = time
        time = time + increment;
        
        %% Controller
        w = x_kp(1:3);
        q = x_kp(4:7);
         
        %% Detumbling
%         wc_B = [0;0;0];
%         q_match = [0;0;0;1];
        
        %% Pointing
        wc_B = [0;0;0];
        q_match = [1;0;0;0];
        
        %% Spin Up
%         Q_BN =SpinCalc('QtoDCM',q',1e-6,0);
%         wc_B = Q_BN*wc_N;
%         wc_B = wc_N;
%         phi = phi+increment*rate;
%         q_match = [wc_B/norm(wc_B)*sin(phi/2); cos(phi/2)]; 
     
        %% Assemble x_tilda
        w_des = w-wc_B;
        q_des = q-q_match;
        x_tilda = [w-wc_B; q-q_match];

        q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
        w1 = w(1); w2 = w(2); w3 = w(3); 
        
        A = Acomputation(Inertia,q1,q2,q3,q4,w1,w2,w3);
        B = Bcomputation(Inertia);
        
        %% Integral addition
%         F = [A zeros(7,3); zeros(3,3) eye(3,3) zeros(3,4)];
%         G = [B; zeros(3,3)];
%         rank(ctrb(F,G))

        %% LQR
         Q = eye(7);
%         Q = [eye(3,3), zeros(3,4) zeros(3,3); 
%              zeros(4,3) eye(4,4) zeros(4,3); 
%              zeros(3,3) zeros(3,4) 10*eye(3,3)];
        R = eye(3,3);
        [K,~,~] = lqr(A,B,Q,R,0);
%         [K,~,~] = lqr(F,G,Q,R,0);
%         torques = -K*[x_tilda; xi(1:3)]; % Integral addition
        torques = -K*x_tilda;             % Proportional alone
        u(counter,:) = torques;
        
%         u_tilda = -K*chi;
%         K_star_inv = inv([A B; eye(7,7) zeros(size(B))]);
%         y_star = [wc_N; q_match];
%         torques = (K_star_inv(2,2)+)
        
        w_diff(counter,:) = w-wc_B;
        q_diff(counter,:) = x_kp(4:6)/sin(acos(x_kp(7)))-[1;0;0];
        
        x_k = updatedSYS(time, time_prev, x_kp, Inertia, torques);
        x_k = normalize(x_k);
        x_system(counter,:) = x_k';
        t(counter) = time;
        
        %% update for next iteration
        xi = xi + x_tilda(4:7)*increment;
        x_kp = x_k;
        counter = counter+1;
    end
    Q_BN = SpinCalc('QtoDCM',q',1e-6,0)
    a = q(1:3)./sin(acos(q(4)))
    
    figure(1)
    subplot(2,1,1)
    plot(t,x_system(:,1),t,x_system(:,2),t,x_system(:,3));
    title('w')
    legend('w_x','w_y','w_z')
    xlabel('Time [s]')
    ylabel('Spin Rate [rad/s]')
    
    subplot(2,1,2)
    plot(t,x_system(:,4),t,x_system(:,5),t,x_system(:,6),t,x_system(:,7));title('q')
    legend('q_1', 'q_2', 'q_3', 'q_4')
    xlabel('Time [s]')
    ylabel('Quaternion Components')
    
    figure(2)
    subplot(3,1,1)
    plot(t,w_diff(:,1),t,w_diff(:,2),t,w_diff(:,3));
    title('w_{difference}')
    legend('w_x','w_y','w_z')
    xlabel('Time [s]')
    ylabel('Spin Rate Error [rad/s]')
    subplot(3,1,2)
    plot(t,q_diff(:,1),'--',t,q_diff(:,2),t,q_diff(:,3),'--');
    title('a_{difference}')
    legend('a_1','a_2','a_3')
    xlabel('Time [s]')
    ylabel('Axis Error')
    subplot(3,1,3)
    title('Torques')
    plot(t,u(:,1),t,u(:,2),t,u(:,3));
    legend('t_1','t_2','t_3')
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    toc
%     animate(x_system(:,4:6),t)
end

function x_normed = normalize(x)
%% Normalize quaternion
q = x(4:7);
x(4:7) = q/(transpose(q)*q);
x_normed = x;
end

function x_k = updatedSYS(t_k, t_km, x_km, Inertia, torques)
%% Propagate dynamics to state array
    dx = SYSdynamics(x_km, Inertia, torques);
    dt = t_k-t_km;
    x_k = x_km + dx*dt;
end

function dx =  SYSdynamics(x_km, Inertia, torques)
GM = 3.986e14; n = length(x_km);
dx = zeros(size(x_km));
%% precompute necessary variables
    w1 = x_km(1);
    w2 = x_km(2);
    w3 = x_km(3);
    q1 = x_km(4);
    q2 = x_km(5);
    q3 = x_km(6);
    q4 = x_km(7);

%% Compute update to state array
    w = [w1; w2; w3];
    q = [q1 q2 q3 q4]';
    qdot = 0.5*[crs(w) w; -w' 0]*q;
    dw = Inertia\(torques - crs(w)*Inertia*w);
    dx(4:7) = qdot;
    dx(1:3) = dw;
    
end