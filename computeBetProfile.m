function [ flag idEnd Arr ] = computeBetProfile( idStart, sdot )
%计算以(idStart,sdot)为起点的正向积分曲线beta-profile
%idStart：正向积分的起点横坐标
%sdot：与idStart对应的线速度
%flag：积分结果标志
%idEnd：正向积分结束点的横坐标（idStart < idEnd）
%Arr：正向积分的线速度，[idStart:1:idEnd]->Arr[1:1:end]

global MVC A B AlpProEnd amount ds;


%beta-profile击穿MVC容忍值：sdot - mvc > tdMVC 表示击穿MVC.单位：m
tdMVC = 0.0006;


Arr = [];
id = idStart;
while(1)
    %击穿s = se
    if id > amount+1
       flag = 1;
       break;
    end
    %击穿MVC
    if sdot > MVC(id) + tdMVC
       flag = 2;
       break;
    end
    %击穿sdot = 0
    if sdot < 0
       flag = 3;
       break;
    end
%     if id >= AlpProEnd.field(1) && id <= AlpProEnd.field(2) && AlpProEnd.Val(id-AlpProEnd.field(1)+1) < sdot
%        flag = 4;
%        break;
%     end
    Arr = [Arr sdot];
    [alp bet] = AlpBet(A(:,id), B(:,id), sdot);
    %调整第一步正向积分的斜率,param:（0,1）
    if id == idStart && (MVC(id+1)-MVC(id))/ds < bet/sdot
       param = 1; 
       bet = alp+((MVC(id+1)-MVC(id))/ds*sdot-alp)*param;
    end
    sdot = sdot + bet/sdot*ds;
    id = id + 1;
end 

idEnd = id-1;

end

