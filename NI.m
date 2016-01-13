clear all
clc
close all

global amount se sdot0 sdote vmax amax R alphal alphar d r ds MVC A B spArr;

amount = 10000;
se = 2*pi ; %m
sdot0 = 0.5;
sdote = 0.1;
vmax = 18;  %18,200
amax = 22;
R = 0.5; 
d = 0.073; 
r = 0.073;
alphal = 2*pi/3;
alphar = -2*pi/3;

ds = se / amount;

%断点阈值
tdDSC = 7;
%击穿阈值：0表示严格的比较，>0表示一定误差的比较
tdBreakDown = 0.001;
% tdTAG = 1;
% tdSGL = 1e-1;
% tdZero = 1e-4;


%开始计时
tic;



%计算MVC曲线，存储在MVC向量
%计算A(s)sdd+B(s)sd^2的A,B系数，存储在A,B向量
A  = [];
B  = [];
MVC= [];
for id = 0 : 1 : amount
   s = ds*id; 
   ab = AbFun(s);
   a = ab(:,1);
   b = ab(:,2);
   A = [A a];
   B = [B b];
   MVC = [MVC Mvc(a, b)];
end


%计算MVC上的alpha和beta值——最小加速度和最大加速度
Alp = [];
Bet = [];
for id = 1 : 1 : amount+1
    [alp bet] = AlpBet(A(:,id), B(:,id), MVC(id));
    Alp = [Alp alp];
    Bet = [Bet bet];
end


%计算MVC上的转换点
%原则：正向和反向积分一步之内不击穿MVC（比较MVC和alpha,beta-profile斜率）
%在满足原则的MVC点中选择断点，切点和奇点
%断点：MVC曲线斜率 > tdDSC
%切点: 以beta正向积分，以alpha反向积分不立即击穿MVC
%奇点：其余点都是奇点
%spArr:存储id:1->amount+1点的类型。0代表普通，1代表断点，2代表切点，3代表奇点
tag = [];
tagSlope = [];
dsc = [];
dscSlope = [];
sgl = [];
sglSlope = [];
spArr = zeros(1,amount+1);
for id = 2 : 1 : amount
    %计算id点的kmvc1(右边斜率)和kmvc2(左边斜率)
    kmvc1 = (MVC(id)-MVC(id-1))/ds;
    kmvc2 = (MVC(id+1)-MVC(id))/ds;
    kalp = Alp(id)/MVC(id);
    kbet = Bet(id)/MVC(id); 
    flag = 0;
    %正向积分不立即击穿MVC
    if kbet <= kmvc2 + tdBreakDown
        flag = flag + 1;
    else
        if kalp <= kmvc2 + tdBreakDown
            flag = flag + 2;
        end
    end
    %反向积分不立即击穿MVC
    if kalp >= kmvc1 - tdBreakDown
       flag = flag + 3;
    else
        if kbet >= kmvc1 - tdBreakDown
            flag = flag + 5;
        end
    end
    
    %flag == 4,5,6,7表示id点正反向积分不立即击穿MVC
    %flag值代表正反向积分的斜率范围
    %4:以bet向前                          ；alp向后
    %5:以(alp, kmvc2)向前一段后以bet向前   ；alp向后
    %6：以bet向前                         ；以(kmvc1,bet)向后一段后以alp向后
    %7：以(alp, kmvc2)向前一段后以bet向前  ；以(kmvc1,bet)向后一段后以alp向后
    if flag == 4 || flag == 5 || flag == 6 || flag == 7
       if (abs(kmvc1) > tdDSC) || (abs(kmvc2) > tdDSC)
           dsc = [dsc id];
           dscSlope = [dscSlope flag];
           spArr(id) = 1;
       else
           if flag == 4 
              tag = [tag id];
              tagSlope = [tagSlope flag];
              spArr(id) = 2;
           else 
              if flag == 5 || flag == 6 || flag == 7
                 sgl = [sgl id];
                 sglSlope = [sglSlope flag];
                 spArr(id) = 3;
              end
           end
       end
    end  
end



%计算MVC转换点的正向和反向积分曲线
%连续奇点是一条反向积分曲线
[flag idEnd Arr] = computeAlpProfile(amount+1, sdote);
alpe.field = [idEnd amount+1];
alpe.val = Arr;
[flag idEnd Arr] = computeBetProfile(1, sdot0);
bet0.field = [1 idEnd];
bet0.val = Arr;
alpProfile = [];
betProfile = [];
for id = 1 : 1 : amount+1
    if spArr(id) ~= 0
        if spArr(id) == 1 || spArr(id) == 2 || (spArr(id)==3 && spArr(id-1)~=3 && spArr(id+1)~=3)
           [flag idEnd Arr] = computeAlpProfile(id, MVC(id));
           alpline.field = [idEnd id];
           alpline.val = Arr;
           [flag idEnd Arr] = computeBetProfile(id, MVC(id));
           betline.field = [id idEnd];
           betline.val = Arr;
           alpProfile = [alpProfile alpline];
           betProfile = [betProfile betline];
        else
            if (spArr(id)==3 && spArr(id-1)~=3 && spArr(id+1)==3)
                tail = -1;
                for i = id : 1 : amount
                    if spArr(i) ~= 3
                       break;
                    end
                    tail = i;
                end
                [flag idEnd Arr] = computeAlpProfile(id, MVC(id));
                alpline.field = [idEnd tail];
                alpline.val = [Arr MVC(id+1:tail)];
                [flag idEnd Arr] = computeBetProfile(tail, MVC(tail));
                betline.field = [tail idEnd];
                betline.val = Arr;
                alpProfile = [alpProfile alpline];
                betProfile = [betProfile betline];
            end
        end
    end
end
alpProfile = [alpProfile alpe];


%计算alpha和beta曲线相交情况
%流程：先确定beta，再确定alpha（与beta最后相交的曲线）
ans = [];
betline = bet0;
start = 1;
szAlpProfile = length(alpProfile);
while (1)
   ia = -1;
   crossnode = -1;
   for i = start : 1 : szAlpProfile
       b1 = betline.field(1);
       b2 = betline.field(2);
       a1 = alpProfile(i).field(1);
       a2 = alpProfile(i).field(2);
       if a1 < b2 && b2 < a2 && b1 <= a1
          %[a1 b2]寻找是否存在交点，如果有记录ia = i，crossnode = 交点横坐标
          cmpRes = betline.val(a1-b1+1:b2-b1+1) - alpProfile(i).val(1 : b2-a1+1);
          for j = 1 : 1 : length(cmpRes)-1
              if cmpRes(j) <= 0 && cmpRes(j+1) >= 0
                 ia = i;
                 crossnode = j+a1-1;
                 break;
              end
          end
       end
   end
   
   if ia == -1
      fprintf('No Ans!');
      break;
   end
   
   ans = [ans betline.val(1:crossnode-b1+1) alpProfile(ia).val(crossnode-alpProfile(ia).field(1)+2:end-1)];
   if ia == szAlpProfile
      fprintf('Get Ans!\n');
      break;
   end
   
   betline = betProfile(ia);
   start = ia+1;
end



%规划时间
timePlan = toc
%计算运动时间
timeTrajectory = 0;
for i = 1 : 1 : amount
    timeTrajectory = timeTrajectory + ds/ans(i);
end
timeTrajectory



figure,plot(ds*[0:1:amount],Bet, '-r', ds*[0:1:amount], Alp,'--g', 'LineWidth', 2);
figure,plot(ds*[0:1:amount], MVC, 'cyan',ds*[1:1:amount], ans,'r', [amount amount]*ds, [0 MVC(amount+1)], '--b', 'LineWidth',2);

% hold on;
% for i = 1 : 1 : length(sgl)
%     plot(ds*(sgl(i)-1), MVC(sgl(i)), 'ro');
% end
% for i = 1 : 1 : length(tag)
%     plot(ds*(tag(i)-1), MVC(tag(i)), 'g*');
% end
% for i = 1 : 1 : length(dsc)
%     plot(ds*(dsc(i)-1), MVC(dsc(i)), 'b+');
% end
% 
% for i = 1 : 1 : length(tag)
%     [flag idEnd Arr] = computeBetProfile(tag(i), MVC(tag(i)));
%     plot(([tag(i):1:idEnd]-1)*ds, Arr, 'r');
%     [flag idEnd Arr] = computeAlpProfile(tag(i), MVC(tag(i)));
%     plot(([idEnd:1:tag(i)]-1)*ds, Arr, 'b');
% end
% for i = 1 : 1 : length(sgl)
%     [flag idEnd Arr] = computeBetProfile(sgl(i), MVC(sgl(i)));
%     plot(([sgl(i):1:idEnd]-1)*ds, Arr, 'r');
%     [flag idEnd Arr] = computeAlpProfile(sgl(i), MVC(sgl(i)));
%     plot(([idEnd:1:sgl(i)]-1)*ds, Arr, 'b');
% end
% 
% [flag idEnd Arr] = computeAlpProfile(amount+1, sdote);
%  plot(([idEnd:1:amount+1]-1)*ds, Arr, 'b');
%  [flag idEnd Arr] = computeBetProfile(1, sdot0);
%  plot(([1:1:idEnd]-1)*ds, Arr, 'r');
% hold off;




