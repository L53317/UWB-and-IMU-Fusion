function [qnb, vn, pos, eth] = insupdate(qnb, vn, pos, wm, vm, ts)  % 捷联惯导数值更新算法
    nn = size(wm,1);  nts = nn*ts;
    [phim, dvbm] = cnscl(wm, vm);  % 圆锥误差/划船误差补偿
    eth = earth(pos, vn);  % 地球相关参数计算
    vn1 = vn + rv2m(-eth.wnin*nts/2)*qmulv(qnb,dvbm) + eth.gcc*nts;  % 速度更新
    vn = (vn+vn1)/2;
    pos = pos + [vn(2)/eth.RMh;vn(1)/eth.clRNh;vn(3)]*nts;  vn = vn1;  % 位置更新
    qnb = qmul(rv2q(-eth.wnin*nts), qmul(qnb, rv2q(phim)));  % 姿态更新
    qnb = qnormlz(qnb);