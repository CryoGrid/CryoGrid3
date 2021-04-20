function fraction = getE_fraction(T, wc, start_reduction)

fraction=double(T>0).*(double(wc>=start_reduction) + double(wc<start_reduction).*0.25.*(1-cos(pi().*wc./start_reduction)).^2);
