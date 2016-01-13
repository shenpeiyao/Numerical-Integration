function [ mat ] = J( betal, yital ,betar, yitar )
%计算机器人运动学模型的J矩阵
%xxxl和xxxr表示机器人偏转角
%mat:J矩阵
global R d r;

mat = [
    -cos(betal)/r  -sin(betal)/r  -R*sin(yital)/r;
    sin(betal)/d   -cos(betal)/d  -R*cos(yital)/d-1;
    -cos(betar)/r  -sin(betar)/r  -R*sin(yitar)/r;
    sin(betar)/d   -cos(betar)/d  -R*cos(yitar)/d-1
    ];

end

