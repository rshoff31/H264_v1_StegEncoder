function Obj = extract_object(F,BB_obj,BB_frame)


[lo,to,wo,ho,ro,bo] = extract_position(BB_obj);

[lf,tf,wf,hf,rf,bf] = extract_position(BB_frame);

if lo < lf
    la = lf;
else
    la = lo;
end

if to < tf
    ta = tf;
else
    ta = to;
end

if ro > rf
    ra = rf;
else
    ra = ro;
end

if bo > bf
    ba = bf;
else
    ba = bo;
end

Obj = zeros(ho,wo);

hoff = la-lo+1;
voff = ta-to+1;
max_width = ra-ro+wo;
max_height = ba-bo+ho;

Obj(voff:max_height,hoff:max_width) = F(ta-tf+1:ba-tf+1,la-lf+1:ra-lf+1);



