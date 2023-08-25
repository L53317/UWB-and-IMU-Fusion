function qpb = qaddphi(qnb, phi)  % 四元数计算值=四元数真实值+失准角误差
    qpb = qmul(rv2q(-phi),qnb);
