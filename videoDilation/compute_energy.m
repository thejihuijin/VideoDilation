% COMPUTE_ENERGY Converts heat maps to energy of frame.
% Frames with higher energy will be slowed down and frames with lower
% energy will be sped up.
% Heat Maps are assumed to be 3D matrices with the same number of "frames".
% Size of individual frames may differ between different heat maps
%
% OF_mags : Optical Flow Magnitudes, assumed to be same size as original
%           video
% saliencyMaps : Saliency of each frame
% tsaliencyMaps : Time-Weighted saliency map of each frame
% method : string determining which heat map to use. Options:
%          - 'OF' - Optical Flow Magnitudes
%          - 'TSAL' - Time-Weighted Saliency Maps
%          - 'MOF' - Saliency Masked Optical Flow
% pool : String determining which pooling function to use. Options:
%          - 'MINK' - Minkowski Pooling, p = 2
%          - 'WP' - Weighted Pooling, p = 1/2
%          - 'FNS' - Five Number Summary
%
% normalized_energy : Energy calculated from heat maps normalized from 0 to
%           1. Of dimension (1, n_frames)

function [normalized_energy] = compute_energy(OF_mags, saliencyMaps, tsaliencyMaps, method, pool)

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Check Parameters
if ~exist('method','var')
    method = 'OF';
end
if ~exist('pool','var')
    pool = 'WP';
end

% Extract Dimensions
[rows, cols, n_frames] = size(OF_mags);

% Compute Energy Maps depending on method
if strcmp(method,'OF')
    % Use Optical Flow Magnitudes
    energy_map = OF_mags;
elseif strcmp(method,'TSAL')
    % Use time weighted saliency maps
    energy_map = tsaliencyMaps;
elseif strcmp(method,'MOF')
    % Use Saliency Masked Optical Flow
    masked_flow_mags = zeros(rows,cols,n_frames);
    for i = 1:n_frames
        % Compute Saliency Mask
        sal_mask = imresize(saliencyMaps(:,:,i),[rows, cols],'bilinear');
        sal_mask = imbinarize(sal_mask);
        masked_flow_mags(:,:,i) = OF_mags(:,:,i).*sal_mask;
    end
    energy_map = masked_flow_mags;
else
    fprintf('Invalid method given.\n')
    fprintf('Options:\n\t- OF - Optical Flow \n');
    fprintf('\t- TSAL - Time Weighted Saliency\n');
    fprintf('\t- MOF - Saliency Masked Optical Flow\n');
    return;
end

% Compute Energy Function depending on pool
if strcmp(pool,'MINK')
    % Minkowski, default p=2
    normalized_energy = minkowski(energy_map);
elseif strcmp(pool,'WP')
    % Weighted Pooling, default p=1/2
    normalized_energy = weight_pool(energy_map);
elseif strcmp(pool,'FNS')
    % Five Num Sum
    normalized_energy = five_num_sum(energy_map);
else
    fprintf('Invalid pooling function given.\n');
    fprintf('Options:\n\t- MINK - Minkowski, p=2\n');
    fprintf('\t- WP - Weighted Pooling, p=1/2\n');
    fprintf('\t- FNS - Five Num Sum\n');
    return;
end

end

