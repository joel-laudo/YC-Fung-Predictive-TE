clear
close all
clc

all_skin_coords_indices_input = readmatrix("P14_INPUT_Skin_Nodes_V1.txt");

all_skin_coords = all_skin_coords_indices_input(:,2:4);
top_layer_skin_coords = all_skin_coords(linspace(1,23789,11895),:);

writematrix(top_layer_skin_coords, "P14_Top_Layer_Skin_Node_Coords.txt")