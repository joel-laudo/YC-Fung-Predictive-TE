function [FixedSkinNodes] = P14_FixSkinNodes_V1(tol, skin_nodes, skin_elements, mesh_pts)

%Generate a list of the skin elements to fix within a certain radius R
skin_nodes(:,5) = zeros(length(skin_nodes(:,1)),1);
skin_elements(:,10) = zeros(length(skin_elements(:,1)),1);
skin_node_coords = skin_nodes(:,2:4);



for i = 1:length(skin_node_coords(:,1))
    skin_node_distances = zeros(length(mesh_pts), 1);
    for j = 1:length(mesh_pts(:,1))
        skin_node_r_vector = skin_node_coords(i, 1:3) - mesh_pts(j, 2:4);
        skin_node_distance = norm(skin_node_r_vector);
        skin_node_distances(j) = skin_node_distance;
    end

    if min(skin_node_distances) > tol
            skin_nodes(i, 5) = 1;
    end
    
end

% Remove zero rows from the skin nodes matrix to only include nodes that
% should be fixed
skin_nodes(all(~skin_nodes(:,5),2), :) = [];

%Generate a list of elements to fix in place

skin_node_list = reshape(skin_nodes(:,1),1,[]);
FixedSkinNodes = skin_node_list;

end