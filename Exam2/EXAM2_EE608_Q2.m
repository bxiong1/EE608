%%
%a)
clf
x = -5:0.05:5;
y = -5:0.05:5;
[X,Y] = meshgrid(x,y);
F = (X-1).^2+2*(Y-2).^2;
v = [0:2:10 10:10:100 100:20:200]
[c,h] = contour(X,Y,F,v,'linewidth',2);
hold on;
axis equal
hf = fill([1 1 -1 -1]*5,[-1 1 1 -1]*5,'w','facealpha',0.8);
feasible = (1-X.^2-Y.^2>0)&(X+Y>0);
F(~feasible)=NaN;
contour(X,Y,F,v,'linewidth',2);

%%
%b) % c)
clc
clear
a = 1;
tol = 0.002;
iter =1;
x(1) = 0.5;
y(1) = 0.5;
h1(1) = compute_h1(x(1),y(1));
h2(1) = compute_h2(x(1),y(1));
val(1) = compute_val(x(1),y(1));

H = compute_hessian(x(1),y(1),a);
d = inv(H);
g = compute_gradient(x(1),y(1),a);
next_val = 0;

while abs(val(iter)-next_val)>tol
    x(iter+1) = x(iter)-d(1,:)*g;
    y(iter+1) = y(iter)-d(2,:)*g;
    h1(iter+1) = compute_h1(x(iter+1),y(iter+1));
    h2(iter+1)  = compute_h2(x(iter+1),y(iter+1));
    val(iter+1) = compute_val(x(iter+1),y(iter+1));
    next_val = val(iter);
    while h1(iter+1)<0||h2(iter+1)<0
        a = a*10
        g = compute_gradient(x(1),y(1),a);
        H = compute_hessian(x(1),y(1),a);
        d = inv(H);
        x(iter+1) = x(iter)-d(1,:)*g;
        y(iter+1) = y(iter)-d(2,:)*g;
        h1(iter+1) = compute_h1(x(iter+1),y(iter+1));
        h2(iter+1)  = compute_h2(x(iter+1),y(iter+1));
        val(iter+1) = compute_val(x(iter+1),y(iter+1));
    	ext_val = val(iter);
    end
    iter = iter+1;
    g = compute_gradient(x(iter),y(iter),a);
    H = compute_hessian(x(iter),y(iter),a);
    d = inv(H);
    a = a*0.5;
end
X = [x(iter),y(iter)]
H1 = compute_h1(X(1),X(2))
H2 = compute_h2(X(1),X(2))
Val = compute_val(X(1),X(2))
for i = 1:iter-1
    X_change = x(i+1)-x(i);
    Y_change = y(i+1)-y(i);
    quiver(x(i),y(i),X_change,Y_change,0)
    hold on
end
function h1 = compute_h1(x,y)
    h1=1-x^2-y^2;
end
function h2 = compute_h2(x,y)
    h2=x+y;
end
function val = compute_val(x,y)
    val = (x-1)^2+2*(y-2)^2;
end
function dx = compute_gradient(x,y,a)
    dx1 = 2*x-2+a*((-2*x/(x^2+y^2-1))-(1/(x+y)));
    dx2 = 4*y-8+a*((-2*y/(x^2+y^2-1))-(1/(x+y)));
    dx = [dx1,dx2]';
end
function H = compute_hessian(x,y,a)
    ddx1 = 2 +((4*a*x^2)/((-x^2-y^2+1)^2)+(2*a/((-x^2-y^2+1)))+a/((x+y)^2));
    dxdy1 = (4*a*x*y/((-x^2-y^2+1)^2))+a/(x+y)^2;
    ddy1 = 4+((4*a*y^2)/((-x^2-y^2+1)^2)+(2*a/((-x^2-y^2+1)))+a/((x+y)^2));
    H = [ddx1, dxdy1;dxdy1,ddy1];
end