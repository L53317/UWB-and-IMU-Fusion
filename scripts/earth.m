function eth = earth(pos, vn)  % 地球相关参数计算
global Re ff wie g0
    ee = sqrt(2*ff-ff^2);  e2 = ee^2;  % 第一偏心率
    eth.sl = sin(pos(1));  eth.cl = cos(pos(1));  eth.tl = eth.sl/eth.cl;%pos(1)为纬度 
    eth.sl2 = eth.sl*eth.sl;  sl4 = eth.sl2*eth.sl2;
    sq = 1-e2*eth.sl2;  sq2 = sqrt(sq);
    eth.RMh = Re*(1-e2)/sq/sq2+pos(3); %子午圈曲率半径(平行纬度)，并添加高度修正
    eth.RNh = Re/sq2+pos(3);           %卯酉圈主曲率半径(平行经度)，并添加高度修正
    eth.clRNh = eth.cl*eth.RNh;
    eth.wnie = wie*[0; eth.cl; eth.sl]; %地球相对惯性系加速度（自转角速度）分量描述
    eth.vn = vn;
    eth.wnen = [-vn(2)/eth.RMh; vn(1)/eth.RNh; vn(1)/eth.RNh*eth.tl];%导航系相对地球系速度
    eth.wnin = eth.wnie + eth.wnen;
    eth.wnien = eth.wnie + eth.wnin;
    gLh = g0*(1+5.27094e-3*eth.sl2+2.32718e-5*sl4)-3.086e-6*pos(3); % grs80重力模型
    eth.gn = [0;0;-gLh];  %导航系当地重力加速度
    eth.gcc = eth.gn - cross(eth.wnien,vn);% 考虑重力/哥氏力/向心力
          %可以使用cross(eth.wnien,vn)代替cros（计算单位时间导航系有害加速度的速度增量）

