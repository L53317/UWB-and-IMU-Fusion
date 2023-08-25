function rv = q2rv(q) % 变换四元数转换为等效旋转矢量
	if q(1)<0,  q = -q;  end
    nmhalf = acos(q(1));  % 等效旋转矢量模值的一半
    if nmhalf>1e-20,   b = 2*nmhalf/sin(nmhalf);
    else             b = 2;    end
    rv = b*q(2:4);   % q = [ cos(|rv|/2); sin(|rv|/2)/|rv|*rv ];