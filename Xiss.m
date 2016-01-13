function [ XISS ] = Xiss( s )
%计算机器人位姿关于s的二阶导
%s:给定路径长度
%XISS：位姿关于s的二阶导

xss = -cos(s/1);
yss = -sin(s/1);
citass = -sin(s);

XISS = [xss yss citass]';

end

