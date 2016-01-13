function [ XIS ] = Xis( s )
%计算机器人运动位姿关于s的一阶导
%s:给定路径的长度
%XIS:位姿关于s的一阶导

xs = -sin(s/1);
ys = cos(s/1);
citas = cos(s);

XIS = [xs ys citas]';

end

