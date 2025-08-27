function proj = proj_quant(x, xq, delta)
% PROJ_QUANT perfoms projection of vector x onto the set of feasible
% solutions for the dequantization problem in the time domain.
%
% Input parameters
%       x       vector of input signal
%       xq      quantized signal
%       delta   quantization step
%
% Pavel Záviška, Brno University of Technology, 2020
% Ondøej Mokrı, Brno University of Technology, 2025

overstep_above = (x - xq) > delta/2;
overstep_below = (xq - x) > delta/2;

proj = x;

proj(overstep_above) = xq(overstep_above) + delta(overstep_above)/2 - eps;
proj(overstep_below) = xq(overstep_below) - delta(overstep_below)/2 + eps;
proj(abs(proj)>1) = 1*sign(proj(abs(proj)>1));

end
