function fraction = getET_fraction(T, wc, start_reduction, zero_reached)

%fraction=double(T>0).*min(1, max(0, (wc-zero_reached)./(start_reduction-zero_reached)));

fraction=double(T>0).*(double(wc>=start_reduction) + double(wc<start_reduction).*0.25.*(1-cos(pi().*wc./start_reduction)).^2);