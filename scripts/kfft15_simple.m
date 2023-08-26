function Ft = kfft15_simple(Cnb, fb)
%简化版的系统矩阵，采用相对位置进行表达
%原型：function Ft = kfft15_simple(eth, Cnb, fb)
%20190429.liu

O33 = zeros(3);I33=eye(3);
%以比力在导航系下的投影为元素构成的反对称阵:fn(Tensor，张量)
fnT = askew(Cnb*fb);
%陀螺仪和加速度计一阶马尔科夫相关时间(3600s)构成的矩阵
betaG=diag([1/3600,1/3600,1/3600]);%陀螺仪
betaA=diag([1/3600,1/3600,1/3600]);%加速度计
Ft = [ O33    O33    O33    -Cnb     O33   %系统矩阵
       fnT    O33    O33     O33     Cnb
       O33    I33    O33     O33     O33
       O33    O33    O33    -betaG   O33
       O33    O33    O33     O33    -betaA];
