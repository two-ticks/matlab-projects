num1 = [10];
den1 =[1 2 10];
num2 = [5];
den2 = [1 5];
% G1(s)= 10/(s^2+2s+10)
% G2(s)= 5/(s+5)
[num, den]=series(num1,den1,num2,den2);
printsys(num,den);
[num, den]=parallel(num1,den1,num2,den2);
printsys(num,den);
[num, den]=feedback(num1,den1,num2,den2);
printsys(num,den);