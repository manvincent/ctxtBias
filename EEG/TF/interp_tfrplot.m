function [TFRinterpRes] = interp_tfrplot(TFRres,interp)
% How much interpolation? Higher value means smoother plot. Test to see
interp_amnt = [];
interp_amnt = interp

% Interpolate  time-frequency power
full_interp_mat = [];
for e = 1:length(TFRres.label)
    mat2d = [];
    mat2d = squeeze(TFRres.powspctrm(e,:,:));
    
    X = [];
    Y = [];
    F = [];
    [X,Y] = meshgrid(1:size(mat2d,2), 1:size(mat2d,1));
    F = scatteredInterpolant(X(:),Y(:),mat2d(:),'linear');
    
    U = [];
    V = [];
    
    [U,V] = meshgrid(linspace(1, size(mat2d,2),interp_amnt), linspace(1,size(mat2d,1),interp_amnt));
        
    p_interp_mat = [];
    p_interp_mat = F(U,V);
    
    full_interp_mat = cat(3,full_interp_mat,p_interp_mat);   
end

TFRres.powspctrm = permute(full_interp_mat,[3 1 2]);

new_time_vec = linspace(min(TFRres.time),max(TFRres.time),interp_amnt);
TFRres.time = [];
TFRres.time = new_time_vec

new_freq_vec = linspace(min(TFRres.freq),max(TFRres.freq),interp_amnt);
TFRres.freq = [];
TFRres.freq = new_freq_vec

TFRinterpRes = [];
TFRinterpRes = TFRres
    