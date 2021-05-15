%%
tic
x = [1,1,1];
epsilon = 1e-6;
alpha = 0.001;
[final_x,iter,norms,val,x_,y_,z_] = run_iterations(x,alpha,epsilon)
toc
final_x
x = final_x(1);
y = final_x(2);
z = final_x(3);
H = [4*y^2+8*z^2+2,8*x*y,16*x*z;8*x*y,4*x^2+2,0;16*x*z,0,8*x^2+2]

is_sym = issymmetric(H)

e_val = eig(H)
isposdef = all(e_val>=0)

iterations = 1:iter;
plot(iterations,val)

figure
plot(iterations,x_)
hold on
plot(iterations,y_)
hold on
plot(iterations,z_)

function dx = compute_gradient(x)
    dx1 = 2*(x(1)+5)+4*x(1)*x(2)^2+8*x(1)*x(3)^2;
    dx2 = 2*(x(2)+8)+4*x(1)^2*x(2);
    dx3 = 2*(x(3)+7)+8*x(1)^2*x(3);
    dx = [dx1,dx2,dx3];
end
function new_x = update_x(x,alpha)
    g=compute_gradient(x);
    new_x = x - alpha*g;
end
function [x,iter,norms,val,x_,y_,z_] = run_iterations(x,alpha,epsilon)
    figure
    iter = 1;
    norms = [];
    val = [];
    x_ = [];
    y_ = [];
    z_ = [];
    x_(end+1) = x(1);
    y_(end+1) = x(2);
    z_(end+1) = x(3);
    gradient = compute_gradient(x);
    f = (x(1)+5)^2+(x(2)+8)^2+(x(3)+7)^2+2*x(1)^2*x(2)^2+4*x(1)^2*x(3)^2;
    val(end+1) = f;
    norms(end+1) = norm(alpha*gradient,2);
    while norms(end)>=epsilon
        x = update_x(x,alpha);
        x_(end+1) = x(1);
        y_(end+1) = x(2);
        z_(end+1) = x(3);
        new_gradient = compute_gradient(x);
        norms(end+1) = norm(new_gradient,2);
        f = (x(1)+5)^2+(x(2)+8)^2+(x(3)+7)^2+2*x(1)^2*x(2)^2+4*x(1)^2*x(3)^2;
        val(end+1) = f;
        iter = iter +1;
    end
end

