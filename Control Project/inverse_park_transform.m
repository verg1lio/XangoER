function [i_alpha, i_beta] = inverse_park_transform(Isd, Isq, theta_e)
    i_alpha = Isd * cos(theta_e) - Isq * sin(theta_e);
    i_beta = Isd * sin(theta_e) + Isq * cos(theta_e);
end