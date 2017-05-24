function Computations
    syms w1 w2 w3 q1 q2 q3 q4  ...
        I1 I2 I3 I4 I5 I6 tau1 tau2 tau3 real;
    Inertia = [I1 I6 I5;
               I6 I2 I4;
               I5 I4 I3];
           
    % w is w~
    
    torques = [tau1 tau2 tau3]';
    w = [w1; w2; w3];
    q = [q1 q2 q3 q4]';
    qdot = 0.5*[crs(w) w; -w' 0]*q;
    dw = Inertia\(torques - crs(w)*Inertia*(w));
    dx = [dw; qdot]; 

    I1 = 100; I2 = I1; I3 = I1;
    tau1 = 0;
    Inertia = diag([I1,I2,I3]);
    torques = [tau1 tau2 tau3]';
    w = [w1; w2; w3];
    q = [q1 q2 q3 q4]';
    
    A = jacobian(dx,[w1 w2 w3 q1 q2 q3 q4])
    
    F = [A zeros(7,7); eye(7) zeros(7,7)]
%     B = jacobian(dx, torques)
%     
%     matlabFunction(A,'file','Acomputation')
%     matlabFunction(B,'file','Bcomputation')
 %% 
%     syms phi real;
%     v_old = sqrt(3)/3*[-1;-1;1]
%     v_new = sqrt(3)/9*[1;1;-5]
    
%     phi = acos(dot(v_old,v_new))
 
%     a = [sqrt(3)/3;sqrt(3)/3;sqrt(3)/3];
%     
%     phi = pi;
%     q = [a*sin(phi/2); cos(phi/2)];
%     
%     q123 = q(1:3);
%     
%     Q = (q(4)^2-q123'*q123)*eye(3,3) + 2*q123*q123' - 2*q(4)*crs(q123)
%     norm([1/3 2/3 2/3])
%     
%     syms D L M ct st real;
%     Q = [ct 0 -st; 0 1 0; st 0 ct];
%     I = [ (D^2*M)/3, -(D*L*M)/4, 0; -(D*L*M)/4, (L^2*M)/3, 0; 0, 0, M*(D^2/4 + L^2/4)];
%     Product = Q*I*transpose(Q)
%    
%     [D,V] = eig(Q);
%     Q*v_old
%     v_new
%     solve(Q*v_old -v_new)


% clc; clear all;
%     syms m1 m2 x1 x2 y1 y2 xd1 xd2 yd1 yd2 g real;
%     M = diag([m1;m1;m2;m2])
%     A = [x1 y1 0 0; (x1-x2) (y1-y2) (x2-x1) (y2-y1); 0 0 0 1]
%     pinv(A)
%     b = [xd1^2+yd1^2; (xd2-xd1)^2+(yd2-yd1)^2; 0]
%     a = [0; g; 0; g]
%     inv(M)
%     simplify(pinv(A*inv(M)))
%     rdd = simplify(a+inv(M)*pinv(A*inv(M))*(b-A*a))
end