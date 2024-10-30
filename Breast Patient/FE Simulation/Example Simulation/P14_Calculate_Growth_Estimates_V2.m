close all
clear
clc

%%%%REWORK THIS THIS
%% Read in information about the skin nodes from the part file including which nodes are in the top layer
P14_Top_Layer_skin_nodes = linspace(1, 23789, 11895);

P14AllSkinNodalCoords = readmatrix("P14_INPUT_Skin_Nodes_V1.txt");
P14_All_Skin_elements = readmatrix("P14_INPUT_Skin_Elements_V1.txt");

%Determine top layer skin nodal coordinates
P14_Top_Layer_skin_nodal_coords = P14AllSkinNodalCoords(P14_Top_Layer_skin_nodes, :);

%Copy all skin elements twice, one to modify throughout the program and one
%to display the result
P14_Top_Layer_Skin_elements_identifier = P14_All_Skin_elements;
P14TopLayerSkinElemsAndNodes = P14_All_Skin_elements;



%% Determine the element set for all the top layer elements and record it for later use extracting growth data from Abaqus


%Loop over all top layer skin nodes

%Replace all instances of top layer skin nodes in the top elements variable with zeros
for i = 1:length(P14_Top_Layer_skin_nodes)
   P14_Top_Layer_Skin_elements_identifier(P14_Top_Layer_Skin_elements_identifier == P14_Top_Layer_skin_nodes(i)) = 0;
end
 

%ensure that the element index column is correct
P14_Top_Layer_Skin_elements_identifier(:,1) = P14_All_Skin_elements(:,1);


%Identify and Remove all element rows which have less than 4 top nodes (4 zeros)
numZerosPerRow = sum(P14_Top_Layer_Skin_elements_identifier == 0, 2);
rowsWithLessThan4Zeros = numZerosPerRow < 4;
P14_Top_Layer_Skin_elements_identifier(rowsWithLessThan4Zeros,:) = [];


P14TopLayerSkinElemsAndNodes(rowsWithLessThan4Zeros,:) = [];



%Write variable that contains the top layer element set

P14TopLayerElemSet = P14TopLayerSkinElemsAndNodes(:,1);


%Remove node columns 6-9 of the top element nodal definition matrix because no top nodes show up here
P14TopLayerSkinElemsAndNodes(:,6:9) = [];


%Iterate through each element in the top layer
P14SkinElemAreas(:,1) = P14TopLayerElemSet;

for i = 1:length(P14TopLayerSkinElemsAndNodes(:,1))
    %Calculate the coordinates of each of the four vertices in the element
    A = P14AllSkinNodalCoords(P14AllSkinNodalCoords(:,1) == P14TopLayerSkinElemsAndNodes(i,2), 2:4);
    B = P14AllSkinNodalCoords(P14AllSkinNodalCoords(:,1) == P14TopLayerSkinElemsAndNodes(i,3), 2:4);
    C = P14AllSkinNodalCoords(P14AllSkinNodalCoords(:,1) == P14TopLayerSkinElemsAndNodes(i,4), 2:4);
    D = P14AllSkinNodalCoords(P14AllSkinNodalCoords(:,1) == P14TopLayerSkinElemsAndNodes(i,5), 2:4);
    
    %Calculate vectors to compute areas
    AB = A-B;
    AD = A-D;
    BC = B-C;
    BD = B-D;
    %Calculate the area of each quadrilateral and store it in an array
    tri_1_area = 0.5*norm(cross(AB,AD));
    tri_2_area = 0.5*norm(cross(BC,BD));
    quad_area = tri_1_area + tri_2_area;
    
    P14SkinElemAreas(i,2) = quad_area;
end



%% Calculate the area of the top surface of each element in the top skin layer

%Read in the growth for the top skin layer from all simulation results
%files

%Define list of nonconverged simulations
nonconverged_list = [7, 23, 56, 72, 80];
for i = 1:1
    if ismember(i, nonconverged_list)
        i = i;
    else
        Elem_growth_data(:,i) = readmatrix(['P14_thg_elem_data_Sim_', num2str(i) , '_MAP.txt']);
    end
end

%Calculate initial surface area
Initial_area = sum(P14SkinElemAreas(:,2));

%Calculate final surface area for all simulations
for i = 1:1
    if ismember(i, nonconverged_list)
        i = i;
    else
        FinalSimGrowthEstimate(i) = dot(transpose(P14SkinElemAreas(:,2)), Elem_growth_data(:,i)) - Initial_area;
    end
end
FinalSimGrowthEstimate = transpose(FinalSimGrowthEstimate);
