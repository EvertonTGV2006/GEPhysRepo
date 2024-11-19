R = [32566, 12486, 10000, 5331, 2490, 1071, 678.1, 387.3];
T = [0, 20, 25, 40, 60, 85, 100, 120];
%
plot(R, T, 'r.');
hold on;

RMAT = zeros(length(R), length(R));
for i = 1 : length(R)
    RMAT(i, 1) = 1;
    for j = 1 : (length(R)-1)
        RMAT(i, j+1) = R(i).^(j);
    end
end
size(RMAT)
RMAT
size(T)
a = RMAT\transpose(T);
y = zeros(30,1);
x = zeros(30,1);
Rvec = zeros(length(R), 1);

a

for k = 1 : 30
    testR = 1000 * k;
    x(k) = testR;
    for i = 1 : length(R)
        Rvec(i) = testR.^(i-1);
    end
    y(k) = dot(Rvec, a);
end

size(y)
size(x)

plot(x, y, 'b-')
hold off;

whos RMAT
