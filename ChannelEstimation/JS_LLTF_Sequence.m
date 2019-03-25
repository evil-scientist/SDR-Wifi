function [lltfLower, lltfUpper] = JS_LLTF_Sequence()
%lltfSequence L-LTF upper and lower subcarrier sequence
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   The sequence used for 20 MHz L-LTF: equation 20-11 in IEEE Std
%   802.11-2012.
 
%   Copyright 2015-2016 The MathWorks, Inc.

%#codegen

    lltfLower = [1; 1; -1; -1; 1; 1; -1; 1; -1; 1; 1; 1; 1; 1; 1; -1; -1; 1; 1; -1; 1; -1; 1; 1; 1; 1;];
    lltfUpper = [1; -1; -1; 1; 1; -1; 1; -1; 1; -1; -1; -1; -1; -1; 1; 1; -1; -1; 1; -1; 1; -1; 1; 1; 1; 1];   
end