% This file is part of script calculating QFI(t) and state coefficients of 
% given state in cavity with displaced mirror.
%
% The function calculates optimal QFI and corresponding state for given 
% time @t, force @f, coupling constant @g, constants @omega0 and @omegaM 
% using @inputState.
% 
% The function uses iterative algorithm for computation of QFI using 
% variational principle from arXiv:1312.1356v1 modified so that it no 
% longer chooses dominant eigenvector but maximizes an operator G below 
% with some constraints using cvxr module. 
%
% Function iterates until obtaining optimal state (newQfi - previousQfi <
% @accuracy) or after @maxSteps.
%
% Return value is a pair: optimal state and qfi.
%
% This script is responsible for displaying and exporting coefficient
% charts.
%
% Author: Marcin Byra, UW
% email: marcin.byra1@gmail.com
% 09/2018

function [state,qfi,steps] = calculateOptimalQFI(inputState)
                           
global chartsVisibility N omegaM omega0 maxSteps accuracy f g t ...
    initialNbar debug

state = inputState; % to avoid changing global variable by accident
qfi = -Inf; % to avoid reaching accuracy after first iteration
nOperator = diag(0:length(inputState)-1); % particle number operator
rho = transpose(kron(state',state)); % density matrix of pure initial state
realT = t*2*pi/omegaM; % current time * unit 

p=0.6;
rho = p*rho + (1-p)*ones(N,N); % Making rho more noisy to minimize error

for step = 1:maxSteps
    
    LambdaRho = lambda_channel(rho, realT);
    LambdaPrimRho = lambdaprim_channel(rho, realT);
    
%     Old way of constructing L
%     L = SLD(LambdaRho, LambdaPrimRho, inputState);

%     Important observation: try/catch because when t = 1, then cvxr raises
%     error stating that its argument is not positive semidefinite.
%     However, when I change t -> t+0.0000001, there is no error, but all
%     charts are sometimes changed and make no sense (mainly in t=1,2,3...)
    try
        cvx_begin sdp quiet
            % Construction of symmetric logarithmic derivative - new way
            
            % let lflat represent L in basis e of NxN hermitian matrices
            variable lflat(N^2)
            e = hermitianBase(N);

            % let first be equivalent of trace(LambdaPrimRho * L)
            first = 0.0; 
            for i = 1:N^2
                first = first + lflat(i)*real(trace(e(:,:,i)*LambdaPrimRho));
            end
            
            % let second be equivalent of trace(LambdaRho * L^2)
            R=zeros(N^2,N^2);
            for p=1:N^2
                for q=1:N^2
                    R(p,q)=trace(e(:,:,q)*e(:,:,p)*LambdaRho);
                end
            end
            R=(R+R')/2;
            second = real(lflat' * R * lflat); 
            
            % now our L (represented as vector l) is result of:
            maximize( 2 * first - second)
        cvx_end
    catch ME
        warning(ME.message)
        fprintf("\n\t\t\tError with optimization of l when t = %.2f ", t);
        qfi = NaN;
        steps = -1;
        return
    end
    
    % restore "quadratic" L 
    L=zeros(N,N);
    for p=1:N^2
        L=L+lflat(p)*e(:,:,p);
    end

    Lsquared = L*L;
    LambdaLsquared = lambda_channel(Lsquared, realT);
    LambdaPrimL = lambdaprim_channel(L, realT);        
    
    LambdaDagLsquared = LambdaLsquared';
    LambdaPrimDagL = LambdaPrimL';
    
    % Construction of maximized operator, general case from arXiv:1312.1356v1
    G = -LambdaDagLsquared + 2*LambdaPrimDagL;

    % maximize using cvx
    cvx_begin sdp quiet
        variable newRho(length(inputState), length(inputState)) semidefinite
        % taking real part because of imaginary artifacts from numerical 
        % inaccuracies(?)
        maximize( real(trace(G * newRho)) ) 
        subject to
            trace(newRho)==1
            trace(newRho * nOperator) == initialNbar
    cvx_end

%     sometimes newRho is sparse, then we have to use eigs() instead of eig
    try
        [Ve,De] = eig(newRho);
        state = Ve(:,N);
    catch ME
        warning(ME.message)
        fprintf("\t\t\tError with eig() when t = %.2f\n", t);
%         TODO run again but deal here with sparse matrix case
%         [Ve,De] = eigs(newRho)
%         assuming eigenvalues are always in increasing order TODO: check
%         state = Ve(:,1);
        qfi = NaN;
        steps = -1;
        return
    end
    newQfi = trace(LambdaRho*L*L);
    
    
    if abs(qfi - newQfi) < accuracy || step == maxSteps

        rho = newRho;
        qfi = newQfi;
        steps = step;
        if (mod(t*4,1) == 0) % every 0.25 time unit
            if debug
                printDebugFoundState(rho, nOperator, Ve, De, state);
            end

            fh = figure(...
                    'Name', sprintf("t = %f", t),...
                    'visible', chartsVisibility);
            bar(abs(state.^2));
%             ylim([0, 1]);
            title(sprintf('t=%f, g=%f, constrained, qfi=%d', t, g, qfi));
            name = sprintf(['\\figures\\coefficients_dim_%d_g_%03d_t_'...
                '%d_constrained_nbar.jpg'], ...
                length(inputState), uint8(g*100), uint8(realT*100));
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

% Based on current rho and t and all constant parameters, constructs state
% after action of quantum channel lambda. The formula is derived in my
% thesis.
function lambda = lambda_channel(rho, t)
global f g omegaM omega0 N initialNbar
lambda = zeros(N, N);
    for a=1:N
        for b=1:N
            lambda(a,b) = rho(a,b) * ...
                exp(-1j*((b-1) - (a-1))*(1-f)*omega0*t) * ...
                exp(1j* (g*(1-2*f)/omegaM)^2 * ((b-1)^2 - (a-1)^2) * (omegaM*t - sin(omegaM*t))^2 ) * ...
                exp(1j* (g*(1-2*f)/omegaM)^2 * ((b-1) - (a-1)) * initialNbar * sin(omegaM*t)) * ...
                exp(...
                    (g*(1-2*f)/omegaM)^2 * ...
                    ( ...
                        ((b-1)-((b-1)-initialNbar)*exp(-1j*omegaM*t)) * ((a-1)-((a-1)-initialNbar)*exp(1j*omegaM*t)) ...
                        - 0.5*((b-1) - ((b-1)-initialNbar)*exp(-1j*omegaM*t)) * ((b-1) - ((b-1)-initialNbar)*exp(1j*omegaM*t)) ...
                        - 0.5*((a-1) - ((a-1)-initialNbar)*exp(-1j*omegaM*t)) * ((a-1) - ((a-1)-initialNbar)*exp(1j*omegaM*t))... 
                    ) ...
                );
        end
    end
end

% Based on current rho and t and all constant parameters, constructs state
% after action of quantum channel lambdaprim, that is, derivative of lambda
% with respect to parameter f (external force). 
% The derivative is calculated copied from Mathematica implementation.
% TODO: implement here using Symbolic Math package.
function lambdaprim = lambdaprim_channel(rho, t)
global f g omegaM omega0 N initialNbar
nbar = initialNbar;
lambdaprim = zeros(N, N);
    for a=1:N
        for b=1:N
            lambdaprim(a,b) = rho(a,b) * ...
                (- omega0*t*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*(a*1i - b*1i) - (g^2*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*(8*f - 4)*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2 - (g^2*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*(8*f - 4)*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2 - (g^2*nbar*exp(-(g^2*(2*f - 1)^2*(- (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(b + exp(-omegaM*t*1i)*(nbar - b + 1) - 1) + (a + exp(omegaM*t*1i)*(nbar - a + 1) - 1)*(a/2 + (exp(-omegaM*t*1i)*(nbar - a + 1))/2 - 1/2) + (b + exp(omegaM*t*1i)*(nbar - b + 1) - 1)*(b/2 + (exp(-omegaM*t*1i)*(nbar - b + 1))/2 - 1/2)))/omegaM^2)*exp(-omega0*t*(f - 1)*(a*1i - b*1i))*exp(-(g^2*nbar*sin(omegaM*t)*(2*f - 1)^2*(a - b)*1i)/omegaM^2)*exp(-(g^2*(2*f - 1)^2*(sin(omegaM*t) - omegaM*t)^2*((a - 1)^2 - (b - 1)^2)*1i)/omegaM^2)*sin(omegaM*t)*(8*f - 4)*(a - b)*1i)/omegaM^2);
        end
    end
end

% Returns set of N^2 ortonormal hermitian matrices of dimension NxN
% e(:, :, k), k=1..N^2 is ortonormal base of NxN matrices
function ei = hermitianBase(D)
    ei = zeros(D,D,D^2);

    for p = 1:D
        ei(p,p,p)=1;
    end
    iter = D+1;
    for p = 2:D
        for q = 1:p-1
            ei(p,q,iter)= 1/sqrt(2);
            ei(q,p,iter)= 1/sqrt(2);
            iter = iter+1;
        end
    end
    for p= 2:D
        for q= 1:p-1
            ei(p,q,iter) = -1i/sqrt(2);
            ei(q,p,iter) = 1i/sqrt(2);
            iter = iter+1;
        end
    end
end

% Deprecated: old way of construction of symmetric logarithmic derivative
% Based on: Matteo G A Paris. Quantum  Estimation  for  Quantum  Technology
function L = SLD(LambdaRho, LambdaPrimRho, inputState)
    [V,D] = eig(LambdaRho);
    L = zeros(length(inputState), length(inputState));
    for k = 1:length(inputState)
        for l = 1:length(inputState)
            if abs(D(k,k) + D(l,l)) > 1e-6
                L = L + 2/(D(k,k) + D(l,l)) * conj(V(:,k)' * ...
                    LambdaPrimRho * V(:,l)) * kron(V(:,k)', V(:,l));
            end    
        end
    end

end


% Writes debug info about found state to stdout.
function printDebugFoundState(rho, nOperator, Ve, De, newState)
    display(rho)
    fprintf('\t\t Tr(rho) = ');
    display(trace(rho));
    fprintf('\t\t Tr(rho * noperator) = ');
    display(trace(rho * nOperator));
    fprintf('\t\t Eigenvectors:');
    display(Ve);
    fprintf('\t\t Corresponding eigenvalues:');
    display(De);
    fprintf(['\t\t |newState><newState|\n\t\t newState is the '...
        'last eigenvector of new rho\n\t\t(works only when last '...
        'eigenvalue is 1 and the rest is close to 0)']);
    kron(newState',newState)
end