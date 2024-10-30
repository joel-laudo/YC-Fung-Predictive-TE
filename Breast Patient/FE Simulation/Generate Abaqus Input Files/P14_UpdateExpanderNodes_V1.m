function [new_node_9999, NewExpNodes] = P14_UpdateExpanderNodes_V1(original_mesh_nodes, original_mesh_elements, expd_pts_and_normals, offset, theta)

%% Notes
%This program generalizes for any 2D mesh and coordinate system except a
%list of the edge nodes must be HANDPICKED and HARDCODED in order around
%the circle
% CHANGE LAMBDA, EDGE NODE INDICES, and POINT CLOUD CSV BETWEEN ITERATIONS


%% Initialize Parameters and Vertex/Element Matrices
%Scaling parameter for the rigid registration
lambda = 1.0; %0.90

%2D mesh points read from a text file taken from the Abaqus input file
mesh_csv = original_mesh_nodes;

%add z = 0 to all points
mesh_csv(:,4) = 0; 

%Pull out the x, y coordinates from the mesh pts file
mesh_pts = mesh_csv(:,2:4);

%Get number of vertices in the top layer
original_number_of_vertices = length(mesh_csv(:,1));



%2D elements read from a text file taken from the Abaqus input file
mesh_elements_csv = original_mesh_elements;

%number of elements in the top mesh layer
top_layer_element_count = length(mesh_elements_csv(:,1)); 

mesh_elements_csv([length(mesh_elements_csv(:,1))+1:2*length(mesh_elements_csv(:,1))], 1) = mesh_elements_csv(:,1) + length(mesh_elements_csv(:,1));

%define the mesh element vertices initially to just include the top layer
%points
mesh_elements_vertex_indices = mesh_elements_csv([1:top_layer_element_count],2:5);

[original_mesh_element_count, vertex_count] = size(mesh_elements_vertex_indices);

%% Arrange Bottom Layer of Elements
bottom_layer_mesh_elements(:,1) = mesh_elements_vertex_indices(:,4) + original_number_of_vertices;
bottom_layer_mesh_elements(:,2) = mesh_elements_vertex_indices(:,3) + original_number_of_vertices;
bottom_layer_mesh_elements(:,3) = mesh_elements_vertex_indices(:,2) + original_number_of_vertices;
bottom_layer_mesh_elements(:,4) = mesh_elements_vertex_indices(:,1) + original_number_of_vertices;


%Define element vertex coordinates by creating an array with the element
%indices and vertex coordinates

mesh_elements_csv([0.5*length(mesh_elements_csv(:,1))+1:length(mesh_elements_csv(:,1))],2:5) = bottom_layer_mesh_elements;

%% Arrange Edge Elements
%HARDCODE edge element node numbers taken directly from Abaqus
edge_nodes(1,:) = [3 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 5 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 4 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 2 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74];
edge_nodes(2,:) = edge_nodes(1,:) + original_number_of_vertices;
edge_elements(:,1) = (1:length(edge_nodes(1,:))) + 2*original_mesh_element_count;

for i = 1:length(edge_elements)
    if i == length(edge_elements)
        edge_elements(i, 2:5) = [edge_nodes(1,i), edge_nodes(1, i) + original_number_of_vertices, edge_nodes(1, 1) + original_number_of_vertices, edge_nodes(1, 1)];
    else
        edge_elements(i,[2:5]) = [edge_nodes(1,i), edge_nodes(1, i) + original_number_of_vertices, edge_nodes(1, i+1) + original_number_of_vertices, edge_nodes(1, i+1)];
    end
end

mesh_elements_csv(edge_elements(:,1), :) = edge_elements;

%% Flip normal vectors for all elements in the mesh so that they point INSIDE the expander
mesh_elements_csv_copy = mesh_elements_csv;
mesh_elements_csv(:,3) = mesh_elements_csv_copy(:,5);
mesh_elements_csv(:,5) = mesh_elements_csv_copy(:,3);

%redefine mesh element vertices to include all elements
mesh_elements_vertex_indices = mesh_elements_csv(:,2:5);

for i = 1:original_mesh_element_count
    for j = 1:vertex_count
        for k = 3*j-1:3*j+1
            mesh_elements_coords_and_indices(i,k) = mesh_csv(mesh_elements_vertex_indices(i,j),k-3*j+3);
        end
    end
end

%Expander coords read from a text file--All points in one line with , as
%delimeter
expander_coords_and_normals_one_line = expd_pts_and_normals;

final_expander_x_s = expander_coords_and_normals_one_line(:,1);
final_expander_y_s = expander_coords_and_normals_one_line(:,2);
final_expander_z_s = expander_coords_and_normals_one_line(:,3);

final_expander_points = [final_expander_x_s, final_expander_y_s, final_expander_z_s];

%pull out the normals from the point cloud
point_cloud_normals = expander_coords_and_normals_one_line(:,4:6);
average_point_cloud_normal = sum(point_cloud_normals)/length(point_cloud_normals);

%Calculate the centroid of the point cloud
centroid_x = sum(final_expander_x_s)/length(final_expander_x_s);
centroid_y = sum(final_expander_y_s)/length(final_expander_y_s);
centroid_z = sum(final_expander_z_s)/length(final_expander_z_s);

centroid = [centroid_x, centroid_y, centroid_z];



%% Calculates the rotation matrix A for the Rigid Registration
% Given initial and final vectors 
A = [0; 0; 1];
B = average_point_cloud_normal;

% Normalize vectors
A = A / norm(A);
B = B / norm(B);

% Calculate rotation axis
R = cross(A, B);

% Calculate rotation angle
cos_angle = dot(A, B);
angle = acos(cos_angle);

% Construct rotation matrix
R_matrix = eye(3) + sin(angle) * skew_symmetric(R) + (1 - cos(angle)) * skew_symmetric(R)^2;

% Construct Scaling Matrix

Scaling_Matrix = [lambda 0 0 ; 0 lambda 0; 0 0 1];
transformed_mesh_points = [];


%% Construct Rotation Matrix

% Construct the rotation matrix
R = [cosd(theta) -sind(theta) 0; sind(theta)  cosd(theta) 0; 0 0 1];

% Rotate the point cloud and apply rigid registration
for i = 1:length(mesh_pts)
    transformed_mesh_points(1:3,i) = R_matrix*Scaling_Matrix*R*transpose(mesh_pts(i,:)) + transpose(centroid);
end



%% Define Transformed mesh vertices and elements
% (just the rigid registration)
transformed_mesh_elements_coords_and_indices(:,1) = mesh_elements_csv(:,1);
transpose_transformed_mesh(:,1) = mesh_csv(:,1);
transpose_transformed_mesh(:,2:4) = transpose(transformed_mesh_points);

for i = 1:original_mesh_element_count
    for j = 1:vertex_count
        for k = 3*j-1:3*j+1
            transformed_mesh_elements_coords_and_indices(i,k) = transpose_transformed_mesh(mesh_elements_vertex_indices(i,j),k-3*j+3);
        end
    end
end


%% Map Expander onto 3D surface
x = [];
%exponential weight for rbf mapping
alpha = 0.4;
[Mesh_rows, Mesh_Cols] = size(transformed_mesh_points);

%Iterate through each 2D mesh pt
for i = 1:Mesh_Cols
    distance_to_travel = [0 0 0];
    r_vectors = [];
    distances = [];

% for each mesh pt, calculate the Cartesian distance to each 3D pt and the
% distance vector and store them
    for j = 1:length(final_expander_x_s) 
        r_vector = transpose(final_expander_points(j,:)) - transformed_mesh_points(:,i);
        distance_between_pts = norm(r_vector); 
        r_vectors(j, :) = r_vector;
        distances(j,1) = distance_between_pts;
    end
    
%Add up the total distance to travel using a weight function for each mesh
%point
    weights = [];

%Calculate weights for each 3D point
    for j = 1:length(final_expander_x_s)
        weights(j) = exp(-(alpha*(distances(j))));
    end

weights = weights/sum(weights);

    for j = 1:length(final_expander_x_s)
        distance_to_travel = distance_to_travel + r_vectors(j,:).*weights(j);
    end
%Record the final distance to travel
    x(1:3, i) = transformed_mesh_points(:, i) + transpose(distance_to_travel);
end

x = transpose(x);
%Offset x towards (-) or away (+) from the skin surface
x= x + offset*average_point_cloud_normal;

final_mesh_x_s = x(:,1);
final_mesh_y_s = x(:,2);
final_mesh_z_s = x(:,3);


% Define bottom layer as 1mm away from the top layer
final_mesh_bottom_layer = x + average_point_cloud_normal;
final_mesh_bottom_x_s = final_mesh_bottom_layer(:,1);
final_mesh_bottom_y_s = final_mesh_bottom_layer(:,2);
final_mesh_bottom_z_s = final_mesh_bottom_layer(:,3);

%% Define 2D to 3D mapped mesh vertices and elements
% Final Expander Coordinates
% Arrange final mesh indices
final_mapped_mesh_elements_coords_and_indices(:,1) = mesh_elements_csv(:,1);
final_mapped_mesh_points(:,1) = mesh_csv(:,1);
final_mapped_mesh_points([length(mesh_csv(:,1))+1:2*length(mesh_csv(:,1))], 1) = mesh_csv(:,1) + length(mesh_csv(:,1));

% Arrange final mesh coordinates
final_mapped_mesh_points([1:length(mesh_csv(:,1))],2:4) = x;
final_mapped_mesh_points([length(mesh_csv(:,1))+1:2*length(mesh_csv)],2:4) = final_mesh_bottom_layer;

for i = 1:length(mesh_elements_csv(:,1))
    for j = 1:vertex_count
        for k = 3*j-1:3*j+1
            final_mapped_mesh_elements_coords_and_indices(i,k) = final_mapped_mesh_points(mesh_elements_vertex_indices(i,j),k-3*j+3);
        end
    end
end




%% Write Output CSV Files and Create Output Variables needed for the Abaqus Input File

%Create CSV files that store the final mapped points and elements for
%Abaqus
NewExpNodes = final_mapped_mesh_points;
new_node_9999 = x(1,:) + 0.5*average_point_cloud_normal;

function S = skew_symmetric(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
end