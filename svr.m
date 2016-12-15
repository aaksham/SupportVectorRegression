
fID = fopen('SVR_dataset.txt','r');
formatSpec = '%f %f';

A = fscanf(fID, formatSpec);
M = reshape(A, [2,50])' ;

h = 0.5; 
C = 4; 
eps = 0.1 ;

s=size(M);
n=s(1,1);
x=M(:,1);
y=M(:,2);

for i=1:n
    for j=1:n
       normSq = (norm(x(i)-x(j)))^2;
       K(i,j)= exp(- (normSq/(2*h*h)) ); 
    end
end

z=zeros(n,1);

cvx_begin
    variable a(n);
    variable a_star(n);
    maximize((-1*0.5*(a-a_star)'*K*(a-a_star))-eps*sum(a+a_star)+(a-a_star)'*y);
    subject to
        a >= 0;
        a_star >= 0;
        a <= C;
        a_star <= C;
cvx_end

ans=cvx_optval;

tps=50;
x_test=linspace(0,1,tps);
y_test=zeros(size(x_test));
for i=1:tps
    kx=zeros(n,1);
    for j=1:n
        normSq = (norm(x(j)-x_test(i)))^2;
        kx(j)= exp(- (normSq/(2*h*h)) );
    end
    y_test(i)=(a-a_star)'*kx;
end

supp_vec_a=a-a_star;
supp_vecs_x=zeros(n,1);
supp_vecs_y=zeros(n,1);
for i=1:n
    if abs(supp_vec_a(i)) > C-exp(-5)
        supp_vecs_x(i)=x(i);
        supp_vecs_y(i)=y(i);
    end
end
        
plot(x_test,y_test);hold on;
scatter(x,y);hold on;
scatter(supp_vecs_x,supp_vecs_y,'g');