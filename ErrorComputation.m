clear
A = []; B = []; C = []; mesh = [];

for n = 1 : 7
    Size = 16*2^(n-1); % set sizes from 16*16 to 1024*1024
    U = FiniteDifferenceMethod(Size);
    A = [A,U(1:2^(n-1):Size)];
    mesh = [mesh, Size];
end
mesh % mesh returns the sizes we simulate

for n = 1 : (size(A,2)-1)
    B = [B,abs(A(:,n+1)-A(:,n))];
end
Error = [NaN,max(B)] % there is no 'error' for the first simulation

for n = 1 : (size(Error,2)-1)
    C = [C,Error(n+1)/Error(n)];
end
OrderOfConvergence = [NaN,log(C)/(log(1/2))]