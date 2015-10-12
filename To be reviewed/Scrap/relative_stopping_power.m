figure;
plot(energy,100*radiative_pwr./total_pwr,'g');
hold on;
plot(energy,100*collision_pwr./total_pwr,'r');
hold on;
plot(energy,100*total_pwr./total_pwr,'black');
hold off;
