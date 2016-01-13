function [ ab ] = AbFun( s )
%计算路劲长s上不等式约束系数向量
%ab:4*2的矩阵，第一列代表a向量，第二列代表b向量
%-vmax<=a*sd<=vmax,sdot>=0
%-amax<=a*sdd+b*sd^2<=amax

global R alphal alphar;

%计算机器人位姿
XI = Xi(s);
x = XI(1);
y = XI(2);
cita = XI(3);

%计算机器人位姿一阶导
XIS = Xis(s);
xs = XIS(1);
ys = XIS(2);
citas = XIS(3);

%计算机器人位姿两阶导
XISS = Xiss(s);
xss = XISS(1);
yss = XISS(2);
citass = XISS(3);

%计算左右轮子的偏转角,xxxl表示左轮,xxxr表示右轮
betal = atan2(-ys-R*citas*cos(alphal+cita), -xs+R*citas*sin(alphal+cita));
yital = betal - alphal - cita;
betar = atan2(-ys-R*citas*cos(alphar+cita), -xs+R*citas*sin(alphar+cita));
yitar = betar - alphar - cita;


%计算J和Jstar矩阵,参考"NI总结"
JMAT = J(betal, yital, betar, yitar);
JSTAR = Jstar(XIS, JMAT, betal, yital, betar, yitar);

%计算a,b向量,参考"NI总结"
ab = [];
ab = [ab JMAT*XIS];
ab = [ab JSTAR*XIS + JMAT*XISS];
end

