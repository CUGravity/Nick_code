function B = Bcomputation(I)
%BCOMPUTATION
%    B = BCOMPUTATION(I1,I2,I3,I4,I5,I6)
I1 = I(1,1);
I2 = I(2,2);
I3 = I(3,3);
I4 = I(2,3);
I5 = I(1,3);
I6 = I(1,2);
%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    12-May-2017 11:58:37

t2 = I4.^2;
t3 = I1.*t2;
t4 = I5.^2;
t5 = I2.*t4;
t6 = I6.^2;
t7 = I3.*t6;
t10 = I1.*I2.*I3;
t11 = I4.*I5.*I6.*2.0;
t8 = t3+t5+t7-t10-t11;
t9 = 1.0./t8;
t12 = I3.*I6;
t13 = t12-I4.*I5;
t14 = t9.*t13;
t15 = I2.*I5;
t16 = t15-I4.*I6;
t17 = t9.*t16;
t18 = I1.*I4;
t19 = t18-I5.*I6;
t20 = t9.*t19;
B = reshape([t9.*(t2-I2.*I3),t14,t17,0.0,0.0,0.0,0.0,t14,t9.*(t4-I1.*I3),t20,0.0,0.0,0.0,0.0,t17,t20,t9.*(t6-I1.*I2),0.0,0.0,0.0,0.0],[7,3]);