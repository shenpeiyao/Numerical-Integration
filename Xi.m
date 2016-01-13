function [ XI ] = Xi( s )
%计算机器人运动位姿xyw
%s:给定路劲上的长度
%XI:对应于s的机器人位姿xyw

x = 1*cos(s/1);
y = 1*sin(s/1);
cita = sin(s);

XI = [x y cita]';

end

