function [ alp bet ] = AlpBet( a, b, mvc )
%计算给定(s,sd)的最小最大加速度
%a,b:AbFun(s)的返回值,大小4*1
%mvc:泛指线速度sd

global amax;

tdZero = 1e-6;

%计算满足加速度约束的最小加速度和最大加速度
alp = -inf;
bet = inf;
for i = 1 : 1 : 4
    if abs(a(i)) > tdZero
       if a(i) > 0
          alp = max(alp, (-amax-b(i)*mvc*mvc)/a(i));
          bet = min(bet, (amax-b(i)*mvc*mvc)/a(i));
       else if a(i) < 0
               alp = max(alp, (amax-b(i)*mvc*mvc)/a(i));
               bet = min(bet, (-amax-b(i)*mvc*mvc)/a(i));
           end
       end
    end
end


end

