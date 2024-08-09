function R = rotation_matrixXYZ(rot_vec)

% Extract angles and convert to rad
alpha = deg2rad(rot_vec(1));
beta = deg2rad(rot_vec(2));
gamma = deg2rad(rot_vec(3));

% Calculate the elements of the rotation matrix
R = [
    cos(gamma)*cos(beta), cos(gamma)*sin(beta)*sin(alpha) - sin(gamma)*cos(alpha), cos(gamma)*sin(beta)*cos(alpha) + sin(gamma)*sin(alpha);
    sin(gamma)*cos(beta), sin(gamma)*sin(beta)*sin(alpha) + cos(gamma)*cos(alpha), sin(gamma)*sin(beta)*cos(alpha) - cos(gamma)*sin(alpha);
    -sin(beta),            cos(beta)*sin(alpha),                                     cos(beta)*cos(alpha)
];

end
