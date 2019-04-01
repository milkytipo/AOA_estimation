function spec = samv(A, Rs, N)

M = size(A,1);
%sigma = norm(x,'fro')^2 /size(A,1)/N;
sigma = norm(Rs,'fro') /M/N;

for k=1:size(A,2) % initialization
     p(k) = real(A(:,k)'*Rs*A(:,k)) / norm(A(:,k))^4;
end

for iter=1:50
    R = A*diag(p)*A' + sigma*eye(M);
    R1 = R^(-1);
    R2 = R^(-2);
    R3 = R1*Rs*R1;
    
    sigma = trace(R2*Rs)/trace(R2);
    
    for k=1:size(A,2)
        p(k) =  real(p(k)*A(:,k)'*R3*A(:,k) / ( A(:,k)'*R1*A(:,k) ));
    end
    
    %iter
    %semilogy(theta,p/max(p),'r-','linewidth',2);grid;xlim([-80 80]);ylim([1e-5 1]);
    %pause(0.2)
end

spec = p';
