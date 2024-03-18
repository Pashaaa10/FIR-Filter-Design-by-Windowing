%% Pasalau Alexandru, 333AB
ng = 1;
ns = 16;

%% FAZA 1
[M,r,beta,L,alfa] = PS_PRJ_1_Faza_1(ng,ns);
% Parametrii ferestrelor
% voi memora valorile din workspace dupa prima rulare
M = 24;
r = 83.1351;
beta = 3.2541;
L = 2;
alfa = 29.4054;

%a - raspunsurile la impuls folosind functia STEM
%figura 1
%primele 4 ferestre (care n-au un parametru specific)
figure()
sgtitle('FIGURA 1 - Fereastra Triunghiulara, Blackman, Hamming si Hanning');
subplot(2, 2, 1);
stem(triang(M));
title('Fereastra Triunghiulara');
subplot(2, 2, 2);
stem(blackman(M));
title('Fereastra Blackman');
subplot(2, 2, 3);
stem(hamming(M));
title('Fereastra Hamming');
subplot(2, 2, 4);
stem(hanning(M));
title('Fereastra Hanning');

%comentarii:
%figura 1 reprezinta cele 4 ferestre (Triunghiulara, Blackman, Hamming si Hanning)
%ele fiind simetrice si similare cu cele din figura 4.3, singura
%diferenta fiind durata ferestrei (Aici e 24) si faptul ca aici
%semnalul este esantionat, nu continual

%figura 2
%3 Cebisev si 3 Kaiser

figure()
sgtitle('FIGURA 2 - Ferestrele Cebisev si Kaiser');
subplot(2,3,1);
stem(chebwin(M,r));
title('Fereastra Cebisev pentru r=',r);
subplot(2,3,2);
stem(chebwin(M,r-5));
title('Fereastra Cebisev pentru r=',r-5);
subplot(2,3,3);
stem(chebwin(M,r+5));
title('Fereastra Cebisev pentru r=',r+5);

subplot(2,3,4);
stem(kaiser(M,beta));
title('Fereastra Kaiser pentru beta=',beta);
subplot(2,3,5);
stem(kaiser(M,beta-1));
title('Fereastra Kaiser pentru beta=',beta-1);
subplot(2,3,6);
stem(kaiser(M,beta+1));
title('Fereastra Kaiser pentru beta=',beta+1);

%comentarii:
%pentru fereastra cebisev, se observa o foarte mica atenuare dar nu in varf
%pentru valorile mai mari ale lui r
%pentru fereastra kaiser, se observa acelasi lucru pentru beta mai mare

%figura 3
%3 Lanczos si 3 Tuckey

figure()
sgtitle('FIGURA 3 - Ferestrele Lanczos si Tuckey');
subplot(2,3,1);
stem(lanczos(M,L));
title('Fereastra Lanczos pentru L=',L);
subplot(2,3,2);
stem(lanczos(M,L-1));
title('Fereastra Lanczos pentru L=',L-1);
subplot(2,3,3);
stem(lanczos(M,L+1));
title('Fereastra Lanczos pentru L=',L+1);

subplot(2,3,4);
stem(tukeywin(M,alfa/100));
title('Fereastra Tuckey pentru alfa=',alfa/100);
subplot(2,3,5);
stem(tukeywin(M,(alfa-15)/100));
title('Fereastra Tuckey pentru alfa=',(alfa-15)/100);
subplot(2,3,6);
stem(tukeywin(M,(alfa+15)/100));
title('Fereastra Tuckey pentru alfa=',(alfa+15)/100);

%comentarii
%fereastra lanczos prezinta un comportament parabolic pentru L=1, dar
%pentru L=2 si L=3, este usor hiperbolic
%fereastra tuckey pare sa aiba marginile din ce in ce mai netede cu cat
%alfa creste (se observa ca pentru alfa = 0.144... arata aproape ca si
%fereastra dreptunghiulara)

%% b
%dreptunghiulara
w_dreptunghi=boxcar(M);
w_dreptunghi = w_dreptunghi/sum(w_dreptunghi);
figure()
sgtitle('FIGURA 4 - Spectrul Ferestrei Dreptunghiulare');
[W_dreptunghi, om_dreptunghi] = freqz(w_dreptunghi,1,5000);
plot(om_dreptunghi/pi, mag2db(abs(W_dreptunghi)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

%cele 4

figure()
sgtitle('FIGURA 5 - Spectrul Ferestrelor Triunghiulare, Blackman, Hamming, Hanning');
w_triunghi = triang(M);
w_triunghi = w_triunghi/sum(w_triunghi);
[W_triunghi, om_triunghi] = freqz(w_triunghi,1,5000);
subplot(2,2,1);
plot(om_triunghi/pi, mag2db(abs(W_triunghi)));
title('Spectrul Ferestrei Triunghiulare');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_blackman = blackman(M);
w_blackman = w_blackman/sum(w_blackman);
[W_blackman, om_blackman] = freqz(w_blackman,1,5000);
subplot(2,2,2);
plot(om_blackman/pi, mag2db(abs(W_blackman)));
title('Spectrul Ferestrei Blackman');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_hamming = hamming(M);
w_hamming = w_hamming/sum(w_hamming);
[W_hamming, om_hamming] = freqz(w_hamming,1,5000);
subplot(2,2,3);
plot(om_hamming/pi, mag2db(abs(W_hamming)));
title('Spectrul Ferestrei Hamming');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_hanning = hanning(M);
w_hanning = w_hanning/sum(w_hanning);
[W_hanning, om_hanning] = freqz(w_hanning,1,5000);
subplot(2,2,4);
plot(om_hanning/pi, mag2db(abs(W_hanning)));
title('Spectrul Ferestrei Hanning');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

%Cebisev si Kaiser

figure()
sgtitle('FIGURA 6 - Spectrele Ferestrelor Cebisev si Kaiser');
w_cebisev = chebwin(M,r);
w_cebisev = w_cebisev/sum(w_cebisev);
[W_cebisev, om_cebisev] = freqz(w_cebisev,1,5000);
subplot(2,3,1);
plot(om_cebisev/pi, mag2db(abs(W_cebisev)));
title('Spectrul Ferestrei Cebisev pentru r=',r);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_cebisev = chebwin(M,r-5);
w_cebisev = w_cebisev/sum(w_cebisev);
[W_cebisev, om_cebisev] = freqz(w_cebisev,1,5000);
subplot(2,3,2);
plot(om_cebisev/pi, mag2db(abs(W_cebisev)));
title('Spectrul Ferestrei Cebisev pentru r=',r-5);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_cebisev = chebwin(M,r+5);
w_cebisev = w_cebisev/sum(w_cebisev);
[W_cebisev, om_cebisev] = freqz(w_cebisev,1,5000);
subplot(2,3,3);
plot(om_cebisev/pi, mag2db(abs(W_cebisev)));
title('Spectrul Ferestrei Cebisev pentru r=',r+5);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_kaiser = kaiser(M,beta);
w_kaiser = w_kaiser/sum(w_kaiser);
[W_kaiser, om_kaiser] = freqz(w_kaiser,1,5000);
subplot(2,3,4);
plot(om_kaiser/pi, mag2db(abs(W_kaiser)));
title('Spectrul Ferestrei Kaiser pentru beta=',beta);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_kaiser = kaiser(M,beta-1);
w_kaiser = w_kaiser/sum(w_kaiser);
[W_kaiser, om_kaiser] = freqz(w_kaiser,1,5000);
subplot(2,3,5);
plot(om_kaiser/pi, mag2db(abs(W_kaiser)));
title('Spectrul Ferestrei Kaiser pentru beta=',beta-1);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_kaiser = kaiser(M,beta+1);
w_kaiser = w_kaiser/sum(w_kaiser);
[W_kaiser, om_kaiser] = freqz(w_kaiser,1,5000);
subplot(2,3,6);
plot(om_kaiser/pi, mag2db(abs(W_kaiser)));
title('Spectrul Ferestrei Kaiser pentru beta=',beta+1);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

%Lanczos si Tuckey

figure()
sgtitle('FIGURA 7 - Spectrele Ferestrelor Lanczos si Tuckey');
w_lanczos = lanczos(M,L);
w_lanczos = w_lanczos/sum(w_lanczos);
[W_lanczos, om_lanczos] = freqz(w_lanczos,1,5000);
subplot(2,3,1);
plot(om_lanczos/pi, mag2db(abs(W_lanczos)));
title('Spectrul Ferestrei Lanczos pentru L=',L);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_lanczos = lanczos(M,L-1);
w_lanczos = w_lanczos/sum(w_lanczos);
[W_lanczos, om_lanczos] = freqz(w_lanczos,1,5000);
subplot(2,3,2);
plot(om_lanczos/pi, mag2db(abs(W_lanczos)));
title('Spectrul Ferestrei Lanczos pentru L=',L-1);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_lanczos = lanczos(M,L+1);
w_lanczos = w_lanczos/sum(w_lanczos);
[W_lanczos, om_lanczos] = freqz(w_lanczos,1,5000);
subplot(2,3,3);
plot(om_lanczos/pi, mag2db(abs(W_lanczos)));
title('Spectrul Ferestrei Lanczos pentru L=',L+1);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_tukey = tukeywin(M,alfa/100);
w_tukey = w_tukey/sum(w_tukey);
[W_tukey, om_tukey] = freqz(w_tukey,1,5000);
subplot(2,3,4);
plot(om_tukey/pi, mag2db(abs(W_tukey)));
title('Spectrul Ferestrei Tuckey pentru alfa=',alfa/100);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_tukey = tukeywin(M,(alfa-15)/100);
w_tukey = w_tukey/sum(w_tukey);
[W_tukey, om_tukey] = freqz(w_tukey,1,5000);
subplot(2,3,5);
plot(om_tukey/pi, mag2db(abs(W_tukey)));
title('Spectrul Ferestrei Tuckey pentru alfa=',(alfa-15)/100);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

w_tukey = tukeywin(M,(alfa+15)/100);
w_tukey = w_tukey/sum(w_tukey);
[W_tukey, om_tukey] = freqz(w_tukey,1,5000);
subplot(2,3,6);
plot(om_tukey/pi, mag2db(abs(W_tukey)));
title('Spectrul Ferestrei Tuckey pentru alfa=',(alfa+15)/100);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');

%comentariile se regasesc in fisierul de prezentare (sub c)

%% FAZA 2
%a
omega_c = PS_PRJ_1_Faza_2a(ng,ns);
omega_c = 0.9238; %memorat din workspace
freq_c = omega_c/pi;
M1 = M-1; %ordinul filtrului

%secventele pondere pentru cele 5 ferestre fara parametru
figure()
sgtitle('FIGURA 8 - Secventele Pondere ale Filtrelor Obtinute Pentru Ferestrele Neparametrice');
dreptunghi_secventa = fir1(M1, freq_c, boxcar(M));
triunghi_secventa = fir1(M1, freq_c, triang(M));
blackman_secventa = fir1(M1, freq_c, blackman(M));
hamming_secventa = fir1(M1, freq_c);
hanning_secventa = fir1(M1, freq_c, hanning(M));

subplot(1, 5, 1);
stem(dreptunghi_secventa);
title('Dreptunghiulara');

subplot(1, 5, 2);
stem(triunghi_secventa);
title('Triunghiulara');

subplot(1, 5, 3);
stem(blackman_secventa);
title('Blackman');

subplot(1, 5, 4);
stem(hamming_secventa);
title('Hamming');

subplot(1, 5, 5);
stem(hanning_secventa);
title('Hanning');

%secventele pondere pentru Cebisev si Kaiser

figure()
sgtitle('FIGURA 9 - Secventele Pondere ale Filtrelor Obtinute Pentru Ferestrele Cebisev si Kaiser');
cebisev_secventa = fir1(M1, freq_c, chebwin(M,r));
cebisev_secventa1 = fir1(M1, freq_c, chebwin(M,r-5));
cebisev_secventa2 = fir1(M1, freq_c, chebwin(M,r+5));
kaiser_secventa = fir1(M1, freq_c, kaiser(M,beta));
kaiser_secventa1 = fir1(M1, freq_c, kaiser(M,beta-1));
kaiser_secventa2 = fir1(M1, freq_c, kaiser(M,beta+1));

subplot(2, 3, 1);
stem(cebisev_secventa);
title('Cebisev cu r=',r);

subplot(2, 3, 2);
stem(cebisev_secventa1);
title('Cebisev cu r=',r-5);

subplot(2, 3, 3);
stem(cebisev_secventa2);
title('Cebisev cu r=',r+5);

subplot(2, 3, 4);
stem(kaiser_secventa);
title('Kaiser cu beta=',beta);

subplot(2, 3, 5);
stem(kaiser_secventa1);
title('Kaiser cu beta=',beta-1);

subplot(2, 3, 6);
stem(kaiser_secventa2);
title('Kaiser cu beta=',beta+1);

%secventele pondere pentru Lanczos si Tuckey

figure()
sgtitle('FIGURA 10 - Secventele Pondere ale Filtrelor Obtinute Pentru Ferestrele Lanczos si Tuckey');
lanczos_secventa = fir1(M1, freq_c, lanczos(M,L));
lanczos_secventa1 = fir1(M1, freq_c, lanczos(M,L-1));
lanczos_secventa2 = fir1(M1, freq_c, lanczos(M,L+1));
tuckey_secventa = fir1(M1, freq_c, tukeywin(M,alfa/100));
tuckey_secventa1 = fir1(M1, freq_c, tukeywin(M,(alfa-15)/100));
tuckey_secventa2 = fir1(M1, freq_c, tukeywin(M,(alfa+15)/100));

subplot(2, 3, 1);
stem(lanczos_secventa);
title('Lanczos cu L=',L);

subplot(2, 3, 2);
stem(lanczos_secventa1);
title('Lanczos cu L=',L-1);

subplot(2, 3, 3);
stem(lanczos_secventa2);
title('Lanczos cu L=',L+1);

subplot(2, 3, 4);
stem(tuckey_secventa);
title('Tuckey cu alfa=',alfa/100);

subplot(2, 3, 5);
stem(tuckey_secventa1);
title('Tuckey cu alfa=',(alfa-15)/100);

subplot(2, 3, 6);
stem(tuckey_secventa2);
title('Tuckey cu alfa=',(alfa+15)/100);

%spectre ferestre fara parametru
figure()
sgtitle('FIGURA 11 - Spectrele si Fazele Filtrelor Obtinute Pentru Ferestrele Neparametrice');
[F_dreptunghi, omega_dreptunghi] = freqz(dreptunghi_secventa,1,5000);
subplot(2,5,1);
plot(omega_dreptunghi/pi, mag2db(abs(F_dreptunghi)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru - Dreptunghiular');

[F_triunghi, omega_triunghi] = freqz(triunghi_secventa,1,5000);
subplot(2,5,2);
plot(omega_triunghi/pi, mag2db(abs(F_triunghi)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru - Triunghiular');

[F_blackman, omega_blackman] = freqz(blackman_secventa,1,5000);
subplot(2,5,3);
plot(omega_blackman/pi, mag2db(abs(F_blackman)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru - Blackman');

[F_hamming, omega_hamming] = freqz(hamming_secventa,1,5000);
subplot(2,5,4);
plot(omega_hamming/pi, mag2db(abs(F_hamming)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru - Hamming');

[F_hanning, omega_hanning] = freqz(hanning_secventa,1,5000);
subplot(2,5,5);
plot(omega_hanning/pi, mag2db(abs(F_hanning)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru - Hanning');

subplot(2,5,6);
plot(omega_dreptunghi/pi, angle(F_dreptunghi));
title('Faza - Dreptunghiular');

subplot(2,5,7);
plot(omega_dreptunghi/pi, angle(F_triunghi));
title('Faza - Triunghiular');

subplot(2,5,8);
plot(omega_dreptunghi/pi, angle(F_blackman));
title('Faza - Blackman');

subplot(2,5,9);
plot(omega_dreptunghi/pi, angle(F_hamming));
title('Faza - Hamming');

subplot(2,5,10);
plot(omega_dreptunghi/pi, angle(F_hanning));
title('Faza - Hanning');

%spectre cebisev si kaiser
figure()
sgtitle('FIGURA 12 - Spectrele si Fazele Filtrelor Obtinute Pentru Ferestrele Cebisev si Kaiser');
[F_cebisev, omega_cebisev] = freqz(cebisev_secventa,1,5000);
subplot(3,4,1);
plot(omega_cebisev/pi, mag2db(abs(F_cebisev)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev pentru r =',r);
subplot(3,4,2);
plot(omega_cebisev/pi, angle(F_cebisev));
title('Faza Cebisev pentru r =',r);

[F_kaiser, omega_kaiser] = freqz(kaiser_secventa,1,5000);
subplot(3,4,3);
plot(omega_kaiser/pi, mag2db(abs(F_kaiser)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser pentru beta =',beta);
subplot(3,4,4);
plot(omega_kaiser/pi, angle(F_kaiser));
title('Faza Kaiser pentru beta =',beta);

[F_cebisev1, omega_cebisev1] = freqz(cebisev_secventa1,1,5000);
subplot(3,4,5);
plot(omega_cebisev1/pi, mag2db(abs(F_cebisev1)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev pentru r =',r-5);
subplot(3,4,6);
plot(omega_cebisev1/pi, angle(F_cebisev1));
title('Faza Cebisev pentru r =',r-5);

[F_kaiser1, omega_kaiser1] = freqz(kaiser_secventa1,1,5000);
subplot(3,4,7);
plot(omega_kaiser1/pi, mag2db(abs(F_kaiser1)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser pentru beta =',beta-1);
subplot(3,4,8);
plot(omega_kaiser1/pi, angle(F_kaiser1));
title('Faza Kaiser pentru beta =',beta-1);

[F_cebisev2, omega_cebisev2] = freqz(cebisev_secventa2,1,5000);
subplot(3,4,9);
plot(omega_cebisev2/pi, mag2db(abs(F_cebisev2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev pentru r =',r+5);
subplot(3,4,10);
plot(omega_cebisev2/pi, angle(F_cebisev2));
title('Faza Cebisev pentru r =',r+5);

[F_kaiser2, omega_kaiser2] = freqz(kaiser_secventa2,1,5000);
subplot(3,4,11);
plot(omega_kaiser2/pi, mag2db(abs(F_kaiser2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser pentru beta =',beta+1);
subplot(3,4,12);
plot(omega_kaiser2/pi, angle(F_kaiser2));
title('Faza Kaiser pentru beta =',beta+1);

%spectre Lanczos si Tuckey
figure()
sgtitle('FIGURA 13 - Spectrele si Fazele Filtrelor Obtinute Pentru Ferestrele Lanczos si Tuckey');
[F_lanczos, omega_lanczos] = freqz(lanczos_secventa,1,5000);
subplot(3,4,1);
plot(omega_lanczos/pi, mag2db(abs(F_lanczos)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Lanczos pentru L =',L);
subplot(3,4,2);
plot(omega_lanczos/pi, angle(F_lanczos));
title('Faza Lanczos pentru L =',L);

[F_tuckey, omega_tuckey] = freqz(tuckey_secventa,1,5000);
subplot(3,4,3);
plot(omega_tuckey/pi, mag2db(abs(F_tuckey)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Tuckey pentru alfa =',alfa/100);
subplot(3,4,4);
plot(omega_tuckey/pi, angle(F_tuckey));
title('Faza Tuckey pentru alfa =',alfa/100);

[F_lanczos1, omega_lanczos1] = freqz(lanczos_secventa1,1,5000);
subplot(3,4,5);
plot(omega_lanczos1/pi, mag2db(abs(F_lanczos1)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Lanczos pentru L =',L-1);
subplot(3,4,6);
plot(omega_lanczos1/pi, angle(F_lanczos1));
title('Faza Lanczos pentru L =',L-1);

[F_tuckey1, omega_tuckey1] = freqz(tuckey_secventa1,1,5000);
subplot(3,4,7);
plot(omega_tuckey1/pi, mag2db(abs(F_tuckey1)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Tuckey pentru alfa =',(alfa-15)/100);
subplot(3,4,8);
plot(omega_tuckey1/pi, angle(F_tuckey1));
title('Faza Tuckey pentru alfa =',(alfa-15)/100);

[F_lanczos2, omega_lanczos2] = freqz(lanczos_secventa2,1,5000);
subplot(3,4,9);
plot(omega_lanczos2/pi, mag2db(abs(F_lanczos2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Lanczos pentru L =',L+1);
subplot(3,4,10);
plot(omega_lanczos2/pi, angle(F_lanczos2));
title('Faza Lanczos pentru L =',L+1);

[F_tuckey2, omega_tuckey2] = freqz(tuckey_secventa2,1,5000);
subplot(3,4,11);
plot(omega_tuckey2/pi, mag2db(abs(F_tuckey2)));
title('Spectru Tuckey pentru alfa =',(alfa+15)/100);
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
subplot(3,4,12);
plot(omega_tuckey2/pi, angle(F_tuckey2));
title('Faza Tuckey pentru alfa =',(alfa+15)/100);

%iau spectrele si le compar

%% clasament
figure()

subplot(4,5,1);
plot(omega_dreptunghi/pi, mag2db(abs(F_dreptunghi)));
title('dreptunghiular');
subplot(4,5,2);
plot(omega_triunghi/pi, mag2db(abs(F_triunghi)));
title('triunghiular');
subplot(4,5,3);
plot(omega_blackman/pi, mag2db(abs(F_blackman)));
title('blackman');
subplot(4,5,4);
plot(omega_hamming/pi, mag2db(abs(F_hamming)));
title('hamming');
subplot(4,5,5);
plot(omega_hanning/pi, mag2db(abs(F_hanning)));
title('hanning');
subplot(4,5,6);
plot(omega_cebisev/pi, mag2db(abs(F_cebisev)));
title('cebisev r');
subplot(4,5,7);
plot(omega_cebisev1/pi, mag2db(abs(F_cebisev1)));
title('cebisev r-5');
subplot(4,5,8);
plot(omega_cebisev2/pi, mag2db(abs(F_cebisev2)));
title('cebisev r+5');
subplot(4,5,9);
plot(omega_kaiser/pi, mag2db(abs(F_kaiser)));
title('kaiser beta');
subplot(4,5,10);
plot(omega_kaiser1/pi, mag2db(abs(F_kaiser1)));
title('kaiser beta-1');
subplot(4,5,11);
plot(omega_kaiser2/pi, mag2db(abs(F_kaiser2)));
title('kaiser beta+1');
subplot(4,5,12);
plot(omega_lanczos/pi, mag2db(abs(F_lanczos)));
title('lanczos L');
subplot(4,5,13);
plot(omega_lanczos1/pi, mag2db(abs(F_lanczos1)));
title('lanczos L-1');
subplot(4,5,14);
plot(omega_lanczos2/pi, mag2db(abs(F_lanczos2)));
title('lanczos L+1');
subplot(4,5,15);
plot(omega_tuckey/pi, mag2db(abs(F_tuckey)));
title('tuckey alfa');
subplot(4,5,16);
plot(omega_tuckey1/pi, mag2db(abs(F_tuckey1)));
title('tuckey alfa-15');
subplot(4,5,17);
plot(omega_tuckey2/pi, mag2db(abs(F_tuckey2)));
title('tuckey alfa+15');

%% 
%am impartit ferestrele pe 3 categorii (in fct de cea mai mica valoare
%afisata in subplot ul anterior) dupa le-am comparat pe categorii

%clasament:
%17 - dreptunghiular
%16 - tuckey -15%
%15 - triunghiular
%14 - kaiser -1
%13 - hamming
%12 - tuckey +15%
%11 - tuckey
%10 - kaiser
%9 - kaiser +1
%8 - lanczos -1
%7 - hanning
%6 - lanczos
%5 - blackman
%4 - lanczos +1
%3 - cebisev -5
%2 - cevisev 
%1 - cebisev +5

%top ferestre:
%1 - cebisev
%2 - blackman
%3 - lanczos
%4 - hanning
%5 - kaiser
%6 - tuckey
%7 - hamming
%8 - triunghiular
%9 - dreptunghiular

%% b
%locul 1
figure()
sgtitle('FIGURA 14 - Spectrele si Fazele Filtrelor - Cebisev cu r+5');
[F_cebisev2, omega_cebisev2] = freqz(cebisev_secventa2,1,5000);
subplot(2,3,1);
plot(omega_cebisev2/pi, mag2db(abs(F_cebisev2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev pentru r = 88.1351 si M = 24');

M_nou1 = M+M/2;
cebisev_secventa2_noua1 = fir1(M_nou1-1, freq_c, chebwin(M_nou1,r+5));
[F_cebisev2_nou1, omega_cebisev2_nou1] = freqz(cebisev_secventa2_noua1,1,5000);
subplot(2,3,2);
plot(omega_cebisev2_nou1/pi, mag2db(abs(F_cebisev2_nou1)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev pentru r = 88.1351 si M = 36');

M_nou2 = 2*M;
cebisev_secventa2_noua2 = fir1(M_nou2-1, freq_c, chebwin(M_nou2,r+5));
[F_cebisev2_nou2, omega_cebisev2_nou2] = freqz(cebisev_secventa2_noua2,1,5000);
subplot(2,3,3);
plot(omega_cebisev2_nou2/pi, mag2db(abs(F_cebisev2_nou2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev pentru r = 88.1351 si M = 48');

subplot(2,3,4);
plot(omega_cebisev2/pi, angle(F_cebisev2));
title('Faza Cebisev pentru r = 88.1351 si M = 24');

subplot(2,3,5);
plot(omega_cebisev2_nou1/pi, angle(F_cebisev2_nou1));
title('Faza Cebisev pentru r = 88.1351 si M = 36');

subplot(2,3,6);
plot(omega_cebisev2_nou2/pi, angle(F_cebisev2_nou2));
title('Faza Cebisev pentru r = 88.1351 si M = 48');

%locul 9
figure()
sgtitle('FIGURA 15 - Spectrele si Fazele Filtrelor - Kaiser cu beta+1');
[F_kaiser2, omega_kaiser2] = freqz(kaiser_secventa2,1,5000);
subplot(2,3,1);
plot(omega_kaiser2/pi, mag2db(abs(F_kaiser2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser pentru beta = 4.2541 si M = 24');

M_nou1 = M+M/2;
kaiser_secventa2_noua1 = fir1(M_nou1-1, freq_c, kaiser(M_nou1,beta+1));
[F_kaiser2_nou1, omega_kaiser2_nou1] = freqz(kaiser_secventa2_noua1,1,5000);
subplot(2,3,2);
plot(omega_kaiser2_nou1/pi, mag2db(abs(F_kaiser2_nou1)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser pentru beta = 4.2541 si M = 36');

M_nou2 = 2*M;
kaiser_secventa2_noua2 = fir1(M_nou2-1, freq_c, kaiser(M_nou2,beta+1));
[F_kaiser2_nou2, omega_kaiser2_nou2] = freqz(kaiser_secventa2_noua2,1,5000);
subplot(2,3,3);
plot(omega_kaiser2_nou2/pi, mag2db(abs(F_kaiser2_nou2)));
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser pentru beta = 4.2541 si M = 48');

subplot(2,3,4);
plot(omega_kaiser2/pi, angle(F_kaiser2));
title('Faza Kaiser pentru beta = 4.2541 si M = 24');

subplot(2,3,5);
plot(omega_kaiser2_nou1/pi, angle(F_kaiser2_nou1));
title('Faza Kaiser pentru beta = 4.2541 si M = 36');

subplot(2,3,6);
plot(omega_kaiser2_nou2/pi, angle(F_kaiser2_nou2));
title('Faza Kaiser pentru beta = 4.2541 si M = 48');

%comentarii
figure();
plot(omega_cebisev2/pi, mag2db(abs(F_cebisev2)));
hold on;
plot(omega_kaiser2_nou2/pi, mag2db(abs(F_kaiser2_nou2)));
title('comparatie pentru observatii');
hold off;
%facut in momentul scrierii documentatiei
%% FAZA 3
%a - functie
%% b

[omega_p,omega_s,Delta_p] = PS_PRJ_1_Faza_3b(ng,ns);
omega_p = 0.9238;
omega_s = 1.0383;
Delta_p = 2.9405;
Delta_s = Delta_p;
M = 24;
sgtitle('FIGURA 16 - Secventele Filtrelor in Ordinea Descrescatoare a Performantei');
subplot(3,3,9);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,7);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,8);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

M_nou1 = 36; % 24 + 24/2 = 36

subplot(3,3,6);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,4);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,5);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

M_nou2 = 72; %24*2 M*3
%filtre imperfecte, a trebuit sa maresc ordinul

subplot(3,3,3);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,1);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,2);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
stem(h);
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);


%% spectrele filtrelor

sgtitle('FIGURA 17 - Spectrele Filtrelor in Ordinea Descrescatoare a Performantei');
subplot(3,3,9);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
hold on;
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,7);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,8);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,6);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,4);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,5);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,3);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,1);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,2);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

%% fazele

sgtitle('FIGURA 18 - Fazele Filtrelor in Ordinea Descrescatoare a Performantei');
subplot(3,3,9);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,7);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,8);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M-1, omega_c / pi, chebwin(M,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=24, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,6);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,4);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,5);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou1-1, omega_c / pi, chebwin(M_nou1,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=36, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,3);
w = 0.25;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,1);
w = 0.5;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

subplot(3,3,2);
w = 0.75;
omega_c = (1-w)*omega_p+w*omega_s;
h = fir1(M_nou2-1, omega_c / pi, chebwin(M_nou2,r+5));
[H, om] = freqz(h,1,5000);
[Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title(['Cebisev, r+5 si M=72, w = ' num2str(w) ', Delta_pr = ' num2str(Delta_pr, '%.2f') ', Delta_sr = ' num2str(Delta_sr, '%.2f')]);

%% FAZA 4
omega_s_mod = 1.2524;
omega_s = omega_s_mod;
%locul 1 GOLD - Kaiser M = 24, beta = 2, omega_c = 1,0937, s = 0,0739
%locul 2 SILVER - Hanning M = 30, omega_c = 1,0937, s = 0,29
%locul 3 BRONZE - Cebisev M = 30, r = 80, omega_c = 1,1, s=0,3533
%locul 4 - Blackman M = 30, omega_c = 1,083, s=0,4198...

sgtitle('FIGURA 19 - Spectre, Faze, Secvente Pondere Pentru Filtrele GOLD, SILVER, BRONZE');
subplot(3,3,1);
M = 24;
beta = 2;
omega_c = 1.0937;
h = fir1(M-1, omega_c / pi, kaiser(M,beta));
[H, om] = freqz(h,1,5000);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser M=24, beta=2, omega_c = 1,0937');

subplot(3,3,2);
M = 24;
omega_c = 1.0937;
h = fir1(M-1, omega_c / pi, hanning(M));
[H, om] = freqz(h,1,5000);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Hanning M=30, omega_c = 1,0937');

subplot(3,3,3);
M = 24;
r = 80;
omega_c = 1.1;
h = fir1(M-1, omega_c / pi, chebwin(M,r));
[H, om] = freqz(h,1,5000);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Cebisev M=30, omega_c = 1,1');

subplot(3,3,4);
M = 24;
beta = 2;
omega_c = 1.0937;
h = fir1(M-1, omega_c / pi, kaiser(M,beta));
[H, om] = freqz(h,1,5000);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title('Faza Kaiser M=24, beta=2, omega_c = 1,0937');

subplot(3,3,5);
M = 24;
omega_c = 1.1;
h = fir1(M-1, omega_c / pi, hanning(M));
[H, om] = freqz(h,1,5000);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title('Faza Hanning M=30, omega_c = 1,0937');

subplot(3,3,6);
M = 24;
r = 80;
omega_c = 1.0837;
h = fir1(M-1, omega_c / pi, chebwin(M,r));
[H, om] = freqz(h,1,5000);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title('Faza Cebisev M=30, r=80, omega_c = 1,1');

subplot(3,3,7);
M = 24;
beta = 2;
omega_c = 1.0937;
h = fir1(M-1, omega_c / pi, kaiser(M,beta));
stem(h);
title('Secventa Pondere Kaiser');

subplot(3,3,8);
M = 24;
omega_c = 1.1;
h = fir1(M-1, omega_c / pi, hanning(M));
stem(h);
title('Secventa Pondere Hanning');

subplot(3,3,9);
M = 24;
r = 80;
omega_c = 1.0837;
h = fir1(M-1, omega_c / pi, chebwin(M,r));
stem(h);
title('Secventa Pondere Cebisev');

save PASALAU_ALEXANDRU_F4 h omega_p omega_c omega_s Delta_p Delta_s;
%% FAZA 5

subplot(1,3,1);
M = 19;
beta = 3;
omega_c = 1.45; %tr sa fie mai mic decat pi/2 (1.57)

omega_p = omega_c-0.34;
omega_s = omega_c+0.35;
Delta_p = 1.9;
Delta_s = 1.34;

%am modificat foarte mult valorile pana am ajuns la filtrul urmator

sgtitle('FIGURA 20 - Spectru, Faza si Secventa Pondere Filtru Optim');
h = fir1(M-1, omega_c / pi, kaiser(M,beta));
[H, om] = freqz(h,1,5000);
plot(om/pi, mag2db(abs(H)));
line([omega_p/pi omega_p/pi], [-200 0], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-200 0], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-200 0], 'Color', 'r');
line([0 omega_p/pi], [mag2db(1+Delta_p/100) mag2db(1+Delta_p/100)], 'Color', 'g');
line([0 omega_p/pi], [mag2db(1-Delta_p/100) mag2db(1-Delta_p/100)], 'Color', 'g');
line([omega_s/pi 1], [mag2db(Delta_s/100) mag2db(Delta_s/100)], 'Color', 'g');
xlabel('Frecventa normalizata (w/pi)');
ylabel('Amplitudine (dB)');
title('Spectru Kaiser M=19, beta=3, omega_c = 1,45');

subplot(1,3,2);
plot(om/pi, angle(H));
line([omega_p/pi omega_p/pi], [-pi pi], 'Color', 'r');
line([omega_c/pi omega_c/pi], [-pi pi], 'Color', 'r');
line([omega_s/pi omega_s/pi], [-pi pi], 'Color', 'r');
title('Faza Kaiser M=19, beta=3, omega_c = 1,45');

subplot(1,3,3);
stem(h);
title('Secventa Pondere Kaiser M=19, beta=3, omega_c = 1,45');
save PASALAU_ALEXANDRU_F#5 h omega_p omega_c omega_s Delta_p Delta_s;