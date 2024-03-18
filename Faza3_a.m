function [Delta_pr, Delta_sr] = Faza3_a(h, omega_p, omega_s, Delta_p, Delta_s)
    k = 0:omega_p/1000:omega_p;
    H = freqz(h, 1, k);
    Delta_pr = max(abs(1-abs(H)));
    k = omega_s:omega_s/1000:pi;
    H = freqz(h,1,k);
    Delta_sr = max(abs(H));
end
