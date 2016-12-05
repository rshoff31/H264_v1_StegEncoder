function [l,t,w,h,r,b] = extract_position(BB)

l = BB(1);
t = BB(2);
w = BB(3);
h = BB(4);
r = l+w-1;
b = t+h-1;