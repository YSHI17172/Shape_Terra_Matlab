function M = HalfEdge(coord, tri)
%HalfEdge Create an efficient half edge type mesh data structure
%   HalfEdge(COORD,TRI) creates a new mesh data structure in half edge
%   format which is very efficient for finding neighbouring points
%   Inputs are XYZ point coordinates COORD and mesh vertices TRI
%   Output is structure M with half edge data structure fields
%   M.vertex
%   M.vertex_incoming_halfedge
%   M.vertex_outgoing_halfedge
%   M.halfedge
%   M.halfedge_next_halfedge
%   M.halfedge_prev_halfedge
%   M.halfedge_opposite_halfedge
%   M.halfedge_parent_edge
%   M.halfedge_parent_facet
%   M.edge
%   M.edge_child_halfedge
%   M.facet
%   M.facet_child_halfedge
%   M.boundary_halfedge_idx

    % Load coord and tri variables
    F = tri;
    nv = size(coord,1);
    nf = size(tri,1);

    % Make a list of non-boundary halfedge + next/prev information
    halfedge_data = reshape(reshape([F(:,1),F(:,2),F(:,2),F(:,3),F(:,3),F(:,1)]', 6*nf, 1), 2, 3*nf)';
    temp_halfedge_nextprev = 1:size(halfedge_data,1)';
    temp_halfedge_nextprev = reshape(temp_halfedge_nextprev, 3, length(temp_halfedge_nextprev)/3);
    halfedge_next_halfedge = [temp_halfedge_nextprev(2,:);temp_halfedge_nextprev(3,:);temp_halfedge_nextprev(1,:)];
    halfedge_next_halfedge = halfedge_next_halfedge(:);
    halfedge_prev_halfedge = [temp_halfedge_nextprev(3,:);temp_halfedge_nextprev(1,:);temp_halfedge_nextprev(2,:)];
    halfedge_prev_halfedge = halfedge_prev_halfedge(:);
    facet_child_halfedge = (1:nf)'*3-2;
    clear temp_halfedge_nextprev

    % Link them to their parent facet
    halfedge_parent_facet = reshape(repmat(1:nf, 3, 1),3*nf,1);

    % Make a list of edges. Keep track on the connection to each halfedge.
    halfedge_order = halfedge_data(:,1) < halfedge_data(:,2);
    temp_edge = [halfedge_order.*halfedge_data(:,1) + (~halfedge_order).*halfedge_data(:,2),...
                 halfedge_order.*halfedge_data(:,2) + (~halfedge_order).*halfedge_data(:,1)];
    temp_edge_child_halfedge = (1:size(temp_edge,1))';
    [edge, ia, ic] = unique(temp_edge,'rows');
    [sorted_ic, iic] = sort(ic);
    halfedge_parent_edge = ic;
    clear halfedge_order temp_edge ia ic

    ne = size(edge,1);

    % Link halfedges to their parent edges. While doing so, create boundary
    % halfedges if there is only one halfedge correspond to a given edge.
    edge_child_halfedge(sorted_ic) = temp_edge_child_halfedge(iic);
    boundary_halfedge = 0;
    i=2;
    while (1)
        if sorted_ic(i) ~= i/2
            sorted_ic = [sorted_ic(1:(i-1));i/2;sorted_ic(i:length(sorted_ic))];
            halfedge_data = [halfedge_data; halfedge_data(edge_child_halfedge(i/2),2) halfedge_data(edge_child_halfedge(i/2),1)];
            halfedge_parent_facet = [halfedge_parent_facet; 0];
            halfedge_next_halfedge = [halfedge_next_halfedge; 0];
            halfedge_prev_halfedge = [halfedge_prev_halfedge; 0];
            iic = [iic(1:(i-1));length(halfedge_data);iic(i:length(iic))];
            halfedge_parent_edge = [halfedge_parent_edge; i/2];
            if boundary_halfedge == 0
                boundary_halfedge = length(halfedge_data);
            end
        end
        i=i+2;
        if (i>length(sorted_ic))
            break
        end
    end
    if mod(length(sorted_ic),2)
        i = length(sorted_ic)+1;
        sorted_ic = [sorted_ic(1:(i-1));i/2];
        halfedge_data = [halfedge_data; halfedge_data(edge_child_halfedge(i/2),2) halfedge_data(edge_child_halfedge(i/2),1)];
        halfedge_parent_facet = [halfedge_parent_facet; 0];
        halfedge_next_halfedge = [halfedge_next_halfedge; 0];
        halfedge_prev_halfedge = [halfedge_prev_halfedge; 0];
        iic = [iic(1:(i-1));length(halfedge_data)];
        halfedge_parent_edge = [halfedge_parent_edge; i/2];
        if boundary_halfedge == 0
            boundary_halfedge = length(halfedge_data);
        end        
    end
    
    edge_child_halfedge = reshape(iic, 2, ne)';
    clear iic temp_edge_child_halfedge

    % Link opposite halfedges
    halfedge_opposite_halfedge(edge_child_halfedge(:,1)) = edge_child_halfedge(:,2);
    halfedge_opposite_halfedge(edge_child_halfedge(:,2)) = edge_child_halfedge(:,1);
    halfedge_opposite_halfedge = halfedge_opposite_halfedge';

    nh = size(halfedge_data,1);

    % Fill out next/prev relations for boundary halfedges
    if boundary_halfedge ~= 0
        for i=boundary_halfedge:nh
            next = halfedge_opposite_halfedge(halfedge_prev_halfedge(halfedge_opposite_halfedge(i)));
            next_prev = halfedge_prev_halfedge(next);
            while next_prev ~= 0
                next = halfedge_opposite_halfedge(halfedge_prev_halfedge(next));
                next_prev = halfedge_prev_halfedge(next);
            end
            halfedge_next_halfedge(i) = next;
            halfedge_prev_halfedge(next) = i;
        end
    end

    [temp_unique_halfedge vertex_outgoing_halfedge] = unique(halfedge_data(:,1));
    [temp_unique_halfedge vertex_incoming_halfedge] = unique(halfedge_data(:,2));

    clear sorted_ic next next_prev i temp_unique_halfedge
    
    
    M.vertex = coord;
    M.vertex_incoming_halfedge = vertex_incoming_halfedge;
    M.vertex_outgoing_halfedge = vertex_outgoing_halfedge;
    M.halfedge = halfedge_data;
    M.halfedge_next_halfedge = halfedge_next_halfedge;
    M.halfedge_prev_halfedge = halfedge_prev_halfedge;
    M.halfedge_opposite_halfedge = halfedge_opposite_halfedge;
    M.halfedge_parent_edge = halfedge_parent_edge;
    M.halfedge_parent_facet = halfedge_parent_facet;
    M.edge = edge;
    M.edge_child_halfedge = edge_child_halfedge;
    M.facet = tri;
    M.facet_child_halfedge = facet_child_halfedge;
    M.boundary_halfedge_idx = boundary_halfedge;
end