%PART 1
%========================%
%Vectors and variables ending in 't' refer to the top spline and those that
%end in 'b' refer to the bottom spline

%x values of nodes
xt = [0.4	0.51	0.69	0.86	1.04	1.34	1.58	2.15	2.75...
    3.72	4.58	6.4	7.94	10.25	12.17	13.92	15.53	17.25...
    18.89	21.4	24.26	27.13	30.08	33.41	37.38 ];
xb = [0.4	0.48	0.58	0.7	0.83	1.15	1.54	2.35	3.08...
    4.23	5.38	7.5	9.34	12.11	14.55	18.43	22.26	24.93...
    27.57	29.7	31.51	33.17	34.77	36.18	37.38];

%y values of nodes
yt = [1.52	1.74	1.92	2.03	2.15	2.29	2.41	2.62	2.82...
    3.08	3.29	3.61	3.8	4.01	4.05	4.08	4.02	3.93...
    3.8	3.57	3.24	2.86	2.41	1.82	0.98];
yb = [1.52	1.42	1.3	1.21	1.13	1.01	0.91	0.77	0.68...
    0.57	0.51	0.44	0.41	0.46	0.53	0.73	0.93	1.06...
    1.14	1.21	1.22	1.21	1.13	1.05	0.98];

%solve for the constant a_j
at = yt;
ab = yb;

%number of nodes
n = length(xt);
%j=n-1

%Solves for and populates h_j
for i = 1:(n-1);
    ht(i) = xt(i+1) - xt(i);
    hb(i) = xb(i+1) - xb(i);
end

%Solves for u_j which is used to populate the diagonal of matrix [A]
for i = 1:(n-2);
    ut(i) = 2 * (ht(i+1) + ht(i));
    ub(i) = 2 * (hb(i+1) + hb(i));
end

%Initializes the matrices [A] and [B] to be used in [A] * [C] = [B]
At = zeros(n-2,n-2);
Bt = zeros(n-2,1);
Ab = zeros(n-2,n-2);
Bb = zeros(n-2,1);

%Populate [A]
for i = 1:(n-2)
    At(i,i) = ut(i);
    Ab(i,i) = ub(i);
end

for i = 1:(n-3)
    
    At(i+1,i)= ht(i+1);
    At(i,i+1)= ht(i+1);
    Ab(i+1,i)= hb(i+1);
    Ab(i,i+1)= hb(i+1);
end

%---------------------%

%Populate [B]
for i = 2:(n-2)
   Bt(i,1) = 3 / ht(i) * (at(i+1) - at(i)) - 3 / ht(i-1) * ...
       (at(i) - at(i-1));
   Bb(i,1) = 3 / hb(i) * (ab(i+1) - ab(i)) - 3 / hb(i-1) * ...
       (ab(i) - ab(i-1));
end

%Solve the system [A] * [Cdummy] = [B] for [Cdummy]
%[Cdummy] = [A]^-1 * [B]
Ainvt = inv(At);
Cdummyt = Ainvt * Bt; 
Ainvb = inv(Ab);
Cdummyb = Ainvb * Bb;

%Populate a new matrix [c]
%Because it is a natural spline, c_1 and c_n = 0;
ct = zeros(n,1);
cb = zeros(n,1);

for i = 1:(n-2);
    ct(i+1) = Cdummyt(i);
    cb(i+1) = Cdummyb(i);
end

%----------------------%

%Solve for dj's
dt = zeros((n-1),1);
db = zeros((n-1),1);

for i = 1:(n-1)
    dt(i) = (ct(i+1) - ct(i)) / (3 * ht(i));
    db(i) = (cb(i+1) - cb(i)) / (3 * hb(i));
end

%Solve for bj's
bt = zeros((n-1),1);
bb = zeros((n-1),1);

for i = 1:(n-1)
bt(i) = (at(i+1) - at(i)) / ht(i) - ht(i) / 3 * (2 * ct(i) + ct(i+1));  
bb(i) = (ab(i+1) - ab(i)) / hb(i) - hb(i) / 3 * (2 * cb(i) + cb(i+1));  
end

%Generate S(j) functions and plot the S(j) spline
figure(1);
hold on;

for i = 1:(n-1)
    tt = linspace(xt(i), xt(i+1), 1000); %generates points for plotting
    %on the appropriate interval for S(j)
    St = at(i) + bt(i) .* (tt - xt(i)) + ct(i) .* ...
        (tt - xt(i)).^2 + dt(i) .* (tt - xt(i)).^3; %S(j)
   
    plot(tt, St);
    
    tb = linspace(xb(i), xb(i+1), 1000);
    Sb = ab(i) + bb(i) .* (tb - xb(i)) + cb(i) .* ...
        (tb - xb(i)).^2 + db(i) .* (tb - xb(i)).^3;
   
    plot(tb, Sb);
end