function [L, M] = Laplacian(coord, tri)
    i1 = tri(:,1); i2 = tri(:,2); i3 = tri(:,3);
    v1 = coord(i3,:) - coord(i2,:);  v2 = coord(i1,:) - coord(i3,:); v3 = coord(i2,:) - coord(i1,:);

    n  = cross(v1,v2,2); 
    dblA = (sqrt(sum((n').^2)))';

    cot12 = -dot(v1,v2,2)./dblA/2; cot23 = -dot(v2,v3,2)./dblA/2; cot31 = -dot(v3,v1,2)./dblA/2;
    diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;

    i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
    j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
    v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
    
    L = sparse(i,j,v,size(coord,1),size(coord,1));

    i = [i1 i2 i3];
    j = [i1 i2 i3];
    diag_v = dblA/6.;
    v = [diag_v,diag_v,diag_v];
    M = sparse(i,j,v,size(coord,1),size(coord,1));
end