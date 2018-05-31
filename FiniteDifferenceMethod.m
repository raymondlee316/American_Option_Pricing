function U = FiniteDifferenceMethod(Size)
% parameter setting
r = 0.01; sigma = 0.2;  K = 100;  T = 1; 
R = 200;
N = Size; M = Size - 1; dt = T/(Size - 1); dx = R/(Size - 1); 
U = zeros(M-1,1); B = zeros(M-1,1); 
A = zeros(M-1,M-1); % the coefficient matrix with a, b, c
BDR = zeros(N-1,1); % store free boundary of option

% Matrix Assignment
for i = 1 : (M-1)
    U(i) = max((K-i*dx),0); % boundary condition at beginning, i.e. U_0_i
end
B = U; V = U;

% Set matrix A
A(1,1) = 1+r*dt+r*1*dt+sigma^2*1^2*dt; % set b1
A(1,2) = -r*1*dt-1/2*sigma^2*1^2*dt; % set c1
A(M-1,M-1) = 1+r*dt+r*(M-1)*dt+sigma^2*(M-1)^2*dt; % set b_M-1
A(M-1,M-2) = -1/2*sigma^2*(M-1)^2*dt; % set a_M-1
for i = 2 : (M-2)
    A(i,i) = 1+r*dt+r*i*dt+sigma^2*i^2*dt; % from the second row, set the diagonal as bi
    A(i,i-1) = -1/2*sigma^2*i^2*dt; % set ai
    A(i,i+1) = -r*i*dt-1/2*sigma^2*i^2*dt; % set ci
end

% Computation
for n = 1 : (N-1)
    B = U;
    B(1) = U(1) -(-1/2*sigma^2*1^2*dt)*K*exp(-r*(n)*dt);
    U = A\B; % B = A*U, so U = inv(A)*B, that is A\B
    
    % Compare half-step price with early excercise payoff
    for i = 1 : size(U)
        if(U(i) < max((K-i*dx),0))
            U(i) = max((K-i*dx),0);
            BDR(n) = i*dx; % Store free boundary
        end
    end
    
    V = [V,U]; % Use matrix V to store the put value at each nodes
end

temp1 = zeros(1,N); % Use temporary matrix to set boundary conditions
for n = 1 : N
    temp1(n) = K*exp(-r*(n-1)*dt); % Put value when x = 0
end
temp2 = zeros(1,N); % Put value when x = R, i.e. stock price is large enough
V = [temp1;V;temp2];

% % Plot the pragh
% % Uncomment to see graphs at each time
% x = 0:dx:200;
% plot(x,V(:,N),'r')
% str = ['The Computed Approximation of American Put When Size=',num2str(Size)];
% title(str)
% str = ['Put Graph When Size=',num2str(Size)];
% print(gcf,'-dpng',str)
% 
% figure
% BDR = flip(BDR);
% x = 0:dt:(1-dt);
% plot(x,BDR,'r')
% str = ['The Free Boundary When Size= ',num2str(Size)];
% title(str)
% str = ['Free Boundary When Size=',num2str(Size)];
% print(gcf,'-dpng',str)

% Return the option value at beginning, i.e. option price
U = V(:,N);