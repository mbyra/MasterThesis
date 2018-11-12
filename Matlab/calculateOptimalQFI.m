% This function calculates optimal QFI and corresponding state for given 
% time @t, force @f, coupling constant @g, constants @omega0 and @omegaM 
% using @inputState.
% 
% The function uses iterative algorithm for computation of QFI using variational
% principle from arXiv:1312.1356v1 modified so that it no longer chooses
% dominant eigenvector but maximizes an operator G below with some
% constraints using cvxr module. 
%
% Function iterates until obtaining optimal state (newQfi - previousQfi <
% @accuracy) or after @maxSteps.
%
% Return value is a pair: optimal state and qfi.
%
%
% Author: Marcin Byra, UW
% email: marcin.byra1@gmail.com
% 08/2018

function [state,qfi,steps] = calculateOptimalQFI(f, g, t, omegaM, omega0, ... 
                               maxSteps, accuracy, inputState, initialNbar, realT)
state = inputState;
qfi = -Inf;
nOperator = diag(0:length(inputState)-1);
nbar = initialNbar;
N = length(inputState);
rho = transpose(kron(state',state));

% Doszumianie rho
p=0.6;
rho = p*rho + (1-p)*ones(N,N);

[R, e] = llrho_flat(rho);

for step = 1:maxSteps
    
    LambdaRho = lambda_channel(rho, t, nbar, f, g, omegaM, omega0, N);
    LambdaPrimRho = lambdaprim_channel(rho, t, nbar, f, g, omegaM, omega0, N);

    [V,D] = eig(LambdaRho);
    
    % Construction of symmetric logarithmic derivative
%     L = zeros(length(inputState), length(inputState));
%     for k = 1:length(inputState)
%         for l = 1:length(inputState)
%             if abs(D(k,k) + D(l,l)) > 1e-6
%                 L = L + 2/(D(k,k) + D(l,l)) * conj(V(:,k)' * LambdaPrimRho...
%                     * V(:,l)) * kron(V(:,k)', V(:,l));
%             end    
%         end
%     end


    % Construction of symmetric logarithmic derivative - new, better way
    
    % Change LambdaRho(N,N) to be (N^2, N^2) in basis e of NxN hermitian matrices
    [R, e] = llrho_flat(LambdaRho);
    D = length(rho);
%     Rprim = zeros(D^2,1);
%     for p=1:D^2
%         Rprim(p) = trace(e(:,:,p) * LambdaPrimRho);
%     end
%     R
%     Rprim
%     Rprim = (Rprim + Rprim')/2;



    % now we can change LambdaRho * L * L to LambdaRho * xflat * xflat
    % where lflat is optimized vector of length N^2
    try
        cvx_begin sdp quiet
            variable lflat(D^2)
    %         maximize( real(2*trace(LambdaPrimRho * L)) - real(trace(LambdaRho * L * L))  )
    %         maximize( real(2 * lflat' * Rprim) - real(lflat' * R * lflat))
            first = 0.0; % equivalent of trace(LambdaPrimRho * L)
            for i = 1:D^2
                first = first + lflat(i)*real(trace(e(:,:,i)*LambdaPrimRho));
            end

            second = real(lflat' * R * lflat); % equivalent of trace(LambdaRho * L^2)

            maximize( 2 * first - second)
        cvx_end
    catch ME
        warning(ME.message)
        fprintf("\t\t\tError with cvx when t = %.2f when finding max L", t);
        rho = newRho;
        qfi = newQfi;
        steps = -1;
        return
    end
    % restore "quadratic" L
    L=zeros(N,N);
    for p=1:N^2
        L=L+lflat(p)*e(:,:,p);
    end

    Lsquared = L*L;
    LambdaLsquared = lambda_channel(Lsquared, t, nbar, f, g, omegaM, omega0, N);
    LambdaPrimL = lambdaprim_channel(L, t, nbar, f, g, omegaM, omega0, N);        
    
    LambdaDagLsquared = LambdaLsquared';
    LambdaPrimDagL = LambdaPrimL';
    
    % Construction of maximized operator, general case from arXiv:1312.1356v1
    G = -LambdaDagLsquared + 2*LambdaPrimDagL;

    % maximize using cvx
    cvx_begin sdp quiet
        variable newRho(length(inputState), length(inputState)) semidefinite
        maximize( real(trace(G * newRho)) ) % imaginary artifacts from numberical inaccuracies
        subject to
            trace(newRho)==1
            trace(newRho * nOperator) == nbar
    cvx_end

    try %there are some exceptions about unability to eig() sparse matrices, ??
        [Ve,De] = eig(newRho);
        newState = Ve(:,N);
    catch ME
        warning(ME.message)
        fprintf("\t\t\tError with eig() when t = %.2f\nRunning again with quiet off:", t);
        [Ve,De] = eigs(newRho)
        newState = Ve(:,1); % assuming eigenvalues are always in increasing order TODO: check
        
    end
    newQfi = trace(LambdaRho*L*L);
    
    
    if abs(qfi - newQfi) < accuracy || step == maxSteps

        rho = newRho;
        qfi = newQfi;
        steps = step;
        display(realT)
        if (mod(realT*4,1) == 0) % every 0.25 time unit
            fprintf('realT = %f, wiec zapisujemy\n', realT);
            fprintf('\t\tFound state = ');
            display(rho)
            fprintf('\t\t Tr(rho) = ');
            display(trace(rho));
            fprintf('\t\t Tr(rho * noperator) = ');
            display(trace(newRho * nOperator));
            Ve
            De
            kron(newState',newState)
            wspolczynnikiKwadrat = abs(newState.^2)
            
            fh = figure('Name', sprintf("t = %f", realT));
            bar(abs(newState.^2));
            title(sprintf('t=%f, g=%f, constrained, qfi=%d', realT, g, qfi));
            name = sprintf('\\figures\\coefficients_dim_%d_g_%03d_t_%d_constrained_nbar.jpg', length(inputState), uint8(g*100), uint8(realT*100));
            saveas(fh,[pwd name]);
%             close(fh);
        end
        return
        
    else
        rho = newRho;
        qfi = newQfi;
        steps = -1;
    end

end
end

function lambda = lambda_channel(rho, t, nbar, f, g, omegaM, omega0, N)
lambda = zeros(N, N);
    for a=1:N
        for b=1:N
            lambda(a,b) = rho(a,b) * ...
                exp(-1j*((b-1) - (a-1))*(1-f)*omega0*t) * ...
                exp(1j* (g*(1-2*f)/omegaM)^2 * ((b-1)^2 - (a-1)^2) * (omegaM*t - sin(omegaM*t))^2 ) * ...
                exp(1j* (g*(1-2*f)/omegaM)^2 * ((b-1) - (a-1)) * nbar * sin(omegaM*t)) * ...
                exp(...
                    (g*(1-2*f)/omegaM)^2 * ...
                    ( ...
                        ((b-1)-((b-1)-nbar)*exp(-1j*omegaM*t)) * ((a-1)-((a-1)-nbar)*exp(1j*omegaM*t)) ...
                        - 0.5*((b-1) - ((b-1)-nbar)*exp(-1j*omegaM*t)) * ((b-1) - ((b-1)-nbar)*exp(1j*omegaM*t)) ...
                        - 0.5*((a-1) - ((a-1)-nbar)*exp(-1j*omegaM*t)) * ((a-1) - ((a-1)-nbar)*exp(1j*omegaM*t))... 
                    ) ...
                );
        end
    end
end

function lambdaprim = lambdaprim_channel(rho, t, nbar, f, g, omegaM, omega0, N)
lambdaprim = zeros(N, N);
    for a=1:N
        for b=1:N
            lambdaprim(a,b) = rho(a,b) * ...
                (- omega0*t*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*(a*1i - b*1i) - (g^2*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*(8*f - 4)*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2 - (g^2*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*(8*f - 4)*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2 - (g^2*nbar*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*sin(omegaM*t)*(8*f - 4)*(a - b)*1i)/omegaM^2);
        end
    end
end

function [R, e] = llrho_flat(rho)
    D = length(rho);

    % Code adopted from W.Górecki
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

    R=zeros(D^2,D^2);
    for p=1:D^2
        for q=1:D^2
            R(p,q)=trace(e(:,:,q)*e(:,:,p)*rho);
        end
    end
    R=(R+R')/2;
end