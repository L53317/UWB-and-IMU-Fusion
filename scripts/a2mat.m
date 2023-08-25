function Cnb = a2mat(att)       
% 姿态角转换为姿态阵，注意方位角北偏西为正
%（符合右手定则，但地理习惯是北偏东为正）
%且规定姿态角的范围为航向角(-pi,pi]，俯仰角[-0.5pi,0.5pi],横滚角(-pi,pi]
%以东北天为导航系xyz，分别对应以右前上为载体系xyz，以符合右手定则转向为正向，
%欧拉角定义为312即：依次绕z0轴正向，x1轴正向，有y2轴正向，得到另一个坐标系(b系)
%%描述导航系(地理坐标系)转到载体坐标系的顺序及转换关系Cnb;当欧拉角均为0时坐标系重合
%注意，这里有vn=Cnb*vb（已经验证过，n上b下）
    s = sin(att); c = cos(att);
    si = s(1); sj = s(2); sk = s(3);   ci = c(1); cj = c(2); ck = c(3);
    Cnb = [ cj*ck-si*sj*sk, -ci*sk,  sj*ck+si*cj*sk;   
            cj*sk+si*sj*ck,  ci*ck,  sj*sk-si*cj*ck;
           -ci*sj,           si,     ci*cj           ];
