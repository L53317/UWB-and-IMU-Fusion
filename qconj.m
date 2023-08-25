function qout = qconj(qin)  % 共轭四元数
    qout = [qin(1); -qin(2:4)];