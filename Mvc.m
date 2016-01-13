function [ mvc ] = Mvc( a, b )
%计算最大速度限制值
%a,b：AbFun(s)返回的不等式约束系数向量,大小4*1
%mvc：最大速度限制值
global vmax amax;

tdZero = 1e-6;

%mvc的计算需要考虑速度和加速度约束
%-vmax<=a*sd<=vmax,sdot>=0
%-amax<=a*sdd+b*sd^2<=amax
mvc = inf;
for i = 1 : 1 : 4
    if abs( a(i) ) > tdZero
        mvc = min(mvc, abs(vmax/a(i)));
        for j = 1 : 1 : 4
           if i ~= j && abs(a(j)) > tdZero
               ci = amax;
               cj = amax;
               if a(i) > 0
                   ci = -ci;
               end
               if a(j) < 0 
                   cj = -cj;
               end
               
               ai = a(i);
               aj = a(j);
               bi = -b(i);
               bj = -b(j);
               
               if abs(bi*aj-bj*ai) > tdZero
                  sdot2 = (cj*ai-ci*aj) / (bi*aj-bj*ai);
                  if sdot2 > 0 
                     mvc = min(mvc, sqrt(sdot2));
                  end
               end
            end
        end
    else if abs(b(i)) > tdZero
            mvc = min(mvc, sqrt( abs(amax/b(i)) ) ); 
        end
    end
end

