function [ mat ] = Jstar( XIS, JS, betal, yital, betar, yitar )
%计算Jstar矩阵，参考"NI总结"

global r d R;

citas = XIS(3);
Jlw = JS(2,:);
Jrw = JS(4,:);

mat = [
       (citas+Jlw*XIS)*sin(betal)/r   -(citas+Jlw*XIS)*cos(betal)/r   -Jlw*XIS*R*cos(yital)/r;
       (citas+Jlw*XIS)*cos(betal)/d   (citas+Jlw*XIS)*sin(betal)/d   Jlw*XIS*R*sin(yital)/d;
       (citas+Jrw*XIS)*sin(betar)/r   -(citas+Jrw*XIS)*cos(betar)/r   -Jrw*XIS*R*cos(yitar)/r;
       (citas+Jrw*XIS)*cos(betar)/d   (citas+Jrw*XIS)*sin(betar)/d   Jrw*XIS*R*sin(yitar)/d;
        ];

end

