%% Matrix exporting

Str = 'Mesh_MATLAB/';

% Reshape data structure for export

sigma_x = reshape(sigma_x,[],1);
sigma_y = reshape(sigma_y,[],1);
sigma_z = reshape(sigma_z,[],1);

dxh = reshape(dxh,[],1);
dyh = reshape(dyh,[],1);
dzh = reshape(dzh,[],1);
dxe = reshape(dxe,[],1);
dye = reshape(dye,[],1);
dze = reshape(dze,[],1);

% Matrix export

writematrix(sigma_x,strcat(Str,'sigmax.txt'))
writematrix(sigma_y,strcat(Str,'sigmay.txt'))
writematrix(sigma_z,strcat(Str,'sigmaz.txt'))

writematrix(dxh,strcat(Str,'dxh.txt'));
writematrix(dyh,strcat(Str,'dyh.txt'));
writematrix(dzh,strcat(Str,'dzh.txt'));

writematrix(dxe,strcat(Str,'dxe.txt'));
writematrix(dye,strcat(Str,'dye.txt'));
writematrix(dze,strcat(Str,'dze.txt'));