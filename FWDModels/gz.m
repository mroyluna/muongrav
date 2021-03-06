function out1 = gz(x1,x2,y1,y2,z1,z2)
%GZ
%    OUT1 = GZ(X1,X2,Y1,Y2,Z1,Z2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    23-Aug-2016 01:30:07

t2 = x1.^2;
t3 = y1.^2;
t4 = z1.^2;
t5 = t2+t3+t4;
t6 = sqrt(t5);
t7 = z2.^2;
t8 = t2+t3+t7;
t9 = sqrt(t8);
t10 = x2.^2;
t11 = t3+t4+t10;
t12 = sqrt(t11);
t13 = y2.^2;
t14 = t2+t4+t13;
t15 = sqrt(t14);
t16 = t3+t7+t10;
t17 = sqrt(t16);
t18 = t2+t7+t13;
t19 = sqrt(t18);
t20 = t4+t10+t13;
t21 = sqrt(t20);
t22 = t7+t10+t13;
t23 = sqrt(t22);
t24 = 1.0./z1;
t25 = 1.0./z2;
out1 = x1.*log(t6+y1).*(-6.67e-11)-y1.*log(t6+x1).*6.67e-11+x1.*log(t9+y1).*6.67e-11+y1.*log(t9+x1).*6.67e-11+x2.*log(t12+y1).*6.67e-11+y1.*log(t12+x2).*6.67e-11+x1.*log(t15+y2).*6.67e-11+y2.*log(t15+x1).*6.67e-11-x2.*log(t17+y1).*6.67e-11-y1.*log(t17+x2).*6.67e-11-x1.*log(t19+y2).*6.67e-11-y2.*log(t19+x1).*6.67e-11-x2.*log(t21+y2).*6.67e-11-y2.*log(t21+x2).*6.67e-11+x2.*log(t23+y2).*6.67e-11+y2.*log(t23+x2).*6.67e-11+z1.*atan(1.0./sqrt(t5).*t24.*x1.*y1).*6.67e-11-z2.*atan(1.0./sqrt(t8).*t25.*x1.*y1).*6.67e-11-z1.*atan(1.0./sqrt(t11).*t24.*x2.*y1).*6.67e-11-z1.*atan(1.0./sqrt(t14).*t24.*x1.*y2).*6.67e-11+z2.*atan(1.0./sqrt(t16).*t25.*x2.*y1).*6.67e-11+z2.*atan(1.0./sqrt(t18).*t25.*x1.*y2).*6.67e-11+z1.*atan(1.0./sqrt(t20).*t24.*x2.*y2).*6.67e-11-z2.*atan(1.0./sqrt(t22).*t25.*x2.*y2).*6.67e-11;
