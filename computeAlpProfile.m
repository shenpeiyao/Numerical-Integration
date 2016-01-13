function [ flag idEnd Arr ] = computeAlpProfile( idStart, sdot )
%计算以(idStart,sdot)为起点的反向积分曲线alpha-profile
%idStart：反向积分的起点横坐标
%sdot：与idStart对应的线速度
%flag：积分结果标志
%idEnd：反向积分结束点的横坐标（idEnd < idStart）
%Arr：反向积分的线速度，[idEnd:1:idStart]->Arr[1:1:end]
global MVC A B ds;

%alpha-profile击穿MVC容忍值：sdot - mvc > tdMVC 表示击穿MVC.单位：m
tdMVC = 0.0006;

Arr = [];
id = idStart;

while(1)
   %击穿s = 0轴
   if id < 1
      flag = 1;
      break;
   end
   %击穿MVC
   if sdot > MVC(id)+tdMVC
      flag = 2;
      break;
   end
   %击穿sdot = 0轴
   if sdot < 0 
      flag = 3;
      break;
   end
   Arr = [sdot Arr];
   [alp bet] = AlpBet(A(:,id), B(:,id), sdot);
   %调整第一步向后积分的斜率,param:（0,1）
   if id == idStart && (MVC(id)-MVC(id-1))/ds > alp/sdot
      param = 1;
      alp = bet+((MVC(id)-MVC(id-1))/ds*sdot-bet)*param;
   end
   sdot = sdot - alp/sdot*ds;
   id = id - 1;
end

idEnd = id + 1;

end

