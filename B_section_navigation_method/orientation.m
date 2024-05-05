%a=[0.99425,7.03703,6.7991];
%h=[8.25,-31.8125,-34];

a=[0.28,-0.108,9.936]
h=[18.75,10.31,-44.43]

%a=[0.0,-0.108,9.98]
%h=[0,20,-44.43]
	
alength=norm(a)
hlength=norm(h)

a=a/norm(a)
h=h/norm(h)

angle=acosd(dot(a,h))

z=-a
x=cross(z,h)
x=x/norm(x)
y=cross(x,z)

R=[x' y' -z']
