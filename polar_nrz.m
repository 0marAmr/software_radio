function [encoded_data] = polar_nrz(data,A)
encoded_data = (2*data-1)*A;
end

