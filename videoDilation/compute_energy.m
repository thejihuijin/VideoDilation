function [normalized_energy] = compute_energy(OF_mags, saliencyMaps, tsaliencyMaps, method, pool)
%COMPUTE_ENERGY

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

