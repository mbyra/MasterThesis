clear;
D=5;
P=2;


rho= rand(D,D)+1i*rand(D,D);
rho= rho+rho';
rho= rho^2;
rho= rho/trace(rho);

state5Constr = [-0.17143-0.483256i, -0.184004-0.466691i,...
    0.0000238363-0.0000768779i, 0.341161+0.277375i, -0.525566-0.125976i];
rho = transpose(kron(state5Constr',state5Constr));
A= rand(D,D)+1i*rand(D,D);
A= A+A';
p=0.6;
rho = p*rho + (1-p)*ones(5,5);
rho
A

%we want to optimize:

%cvx_begin sdp quiet
%variable X(D,D) hermitian;
%minimize(trace(X*X*rho));
%subject to
%    trace((X-A)*(X-A))<=1;
%cvx_end

%but CVX cannot do that

%{e(:,:,k)} is orthonormal base of DxD Hermitian matrices, which will be used in calculations
e= zeros(D,D,D^2);


for p= 1:D
    e(p,p,p)=1;
end
iter= D+1;
for p= 2:D
    for q= 1:p-1
        e(p,q,iter)= 1/sqrt(2);
        e(q,p,iter)= 1/sqrt(2);
        iter= iter+1;
    end
end
for p= 2:D
    for q= 1:p-1
        e(p,q,iter)= -1i/sqrt(2);
        e(q,p,iter)= 1i/sqrt(2);
        iter= iter+1;
    end
end
e

R=zeros(D^2,D^2);
for p=1:D^2
    for q=1:D^2
        R(p,q)=trace(e(:,:,q)*e(:,:,p)*rho);
    end
end
R=(R+R')/2;


aflat=zeros(D^2,1);
for p=1:D^2
    aflat(p,1)=trace(A*e(:,:,p));
end

R


cvx_begin sdp quiet
variable xflat(D^2);
minimize(xflat'*R*xflat);
subject to
    (xflat-aflat)'*(xflat-aflat)<=1;
cvx_end

xflat
xflat'

X=zeros(D,D);
for p=1:D^2
    X=X+xflat(p)*e(:,:,p);
end

display('(xflat-aflat)^T*(xflat-aflat)=');
(xflat-aflat)'*(xflat-aflat)
display('trace((X-A)*(X-A))=');
trace((X-A)*(X-A))
display('trace(X*X*rho)=');
trace(X*X*rho)
display('xflat^T*R*xflat=');
xflat'*R*xflat








