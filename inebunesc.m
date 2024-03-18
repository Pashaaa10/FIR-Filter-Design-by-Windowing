ng = 1;
ns = 16;

[omega_p,omega_s,Delta_p] = PS_PRJ_1_Faza_3b(ng,ns);
Delta_s = Delta_p;
%am decis sa folosesc omega_s_mod, avand valoarea 1.2524 (generata aleator
%in momentul rularii comenzii de pe linia 4), pentru a avea aceeasi banda
%de trecere la fiecare rulare.
omega_s_mod = 1.2524;
omega_c = omega_p:0.01:omega_s_mod;
M = 16:1:30;
r = 80:1:100;

%% Cebisev
for i = 1 :length(M)
    for j = 1 : length(omega_c)
        for k = 1:length(r)
            h = fir1(M(i)-1, omega_c(j) / pi, chebwin(M(i),r(k)));
            [Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s_mod, Delta_p, Delta_s);
            sum = Delta_pr + Delta_sr;
            fprintf(" M = %d, r = %d, omega_c = %f, sum = %f \n",M(i), r(k),omega_c(j),sum);
        end
     end
end

%% Lanczos


L = 2:1:7;
for i = 1 :length(M)
    for j = 1 : length(omega_c)
        for k = 1:length(L)
            h = fir1(M(i)-1, omega_c(j) / pi, lanczos(M(i),L(k)));
            [Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s_mod, Delta_p, Delta_s);
            sum = Delta_pr + Delta_sr;
            fprintf(" M = %d, L = %d, omega_c = %f, sum = %f \n",M(i), L(k),omega_c(j),sum);
        end
     end
end
%crapa codul
%% Blackman

for i = 1 :length(M)
    for j = 1 : length(omega_c)
        h = fir1(M(i)-1, omega_c(j) / pi, blackman(M(i)));
        [Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s_mod, Delta_p, Delta_s);
        sum = Delta_pr + Delta_sr;
        fprintf(" M = %d, omega_c = %f, sum = %f \n",M(i),omega_c(j),sum);
    end
end

%% Hanning

for i = 1 :length(M)
    for j = 1 : length(omega_c)
        h = fir1(M(i)-1, omega_c(j) / pi, hanning(M(i)));
        [Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s_mod, Delta_p, Delta_s);
        sum = Delta_pr + Delta_sr;
        fprintf(" M = %d, omega_c = %f, sum = %f \n",M(i),omega_c(j),sum);
    end
end

%% Kaiser

beta = 2:1:10;
for i = 1 :length(M)
    for j = 1 : length(omega_c)
        for k = 1:length(beta)
            h = fir1(M(i)-1, omega_c(j) / pi, kaiser(M(i),beta(k)));
            [Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s_mod, Delta_p, Delta_s);
            sum = Delta_pr + Delta_sr;
            fprintf(" M = %d, beta = %d, omega_c = %f, sum = %f \n",M(i), beta(k),omega_c(j),sum);
        end
     end
end

%locul 1 - Kaiser M = 24, beta = 2, omega_c = 1,0837, s = 0,1629
%locul 2 - Kaiser M = 30, beta = 2, omega_c = 1,0937, s = 0,0739
%locul 3 - Hanning M = 30, omega_c = 1,0937, s = 0,29
%locul 4 - Cebisev M = 30, r = 80, omega_c = 1,1, s=0,3533
%locul 5 - Blackman M = 30, omega_c = 1,083, s=0,4198