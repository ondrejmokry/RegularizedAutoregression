function [PSM, PSMt, ODG, PSM_inst] = audioqual_silent(RefSig, TestSig, fs) %#ok<*INUSD>

[~, PSM, PSMt, ODG, PSM_inst] = evalc('audioqual(RefSig, TestSig, fs)');

end

