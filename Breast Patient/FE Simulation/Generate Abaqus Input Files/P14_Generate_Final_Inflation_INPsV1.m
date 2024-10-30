clc
clear
close all;

%% Load Base INP file
fileID = fopen('P14_Final_Protocol_Needs_Edit_V1.inp','r');
Incomplete = textscan( fileID, '%s', 'Delimiter','\n' );
Incomplete = Incomplete{1};
fclose(fileID);

%% Mark out the where additional information needs to be added
% AddRowPara=11; %adds k information
AddRowExpNodes=82933; %adds shifted and rotated expander nodes
AddRowNode9999=89510; %adds the new node 9999 after shifting
AddRowSpringStiffness=89645; %adds spring stiffness information
AddRowFixed=90460; %adds fixed skin nodes



%% Load Skin Nodes and Skin Elements
Skin_Nodes = readmatrix("P14_INPUT_Skin_Nodes_V1.txt");
Skin_Elements = readmatrix("P14_INPUT_Skin_Elements_V1.txt");

%% Load original mesh points and mesh elements
original_mesh_nodes = readmatrix("P14_Expander_mesh_nodes_V1.txt");
original_mesh_elements = readmatrix("P14_Expander_mesh_elements_V1.txt");
expd_pts_and_normals = readmatrix("P14_ExpanderSkin_points_with_normals.xyz", "FileType","text");
%% Load experimental design
designs=load('P14_tol_h_theta_Sx_Tcrit_mu_MAP_design.txt');

%% Write Job files,
for JobNum=1:1
    tol = designs(JobNum,2);
    offset = designs(JobNum, 3);
    theta = designs(JobNum, 4);
    Sx = designs(JobNum, 5);
    tcrit = designs(JobNum, 6);
    mu = designs(JobNum, 7);
    fname='FINAL_P14_Job'+string(JobNum)+'_MAP.inp';
    fid = fopen(fname, 'w' );
    
    v = 0.45; %poisson's ratio
    lam = (2*mu*v)/(1-2*v);

    for ii = 1 : 11
        fprintf(fid, '%s\n', Incomplete{ii} );
    end
    %Add mu
    fprintf(fid, "lam = %4.4f\n", lam);
    fprintf(fid, "mu = %4.4f\n", mu);

    for ii = 12 : 16
        fprintf(fid, '%s\n', Incomplete{ii} );
    end

    fprintf(fid, "tcrt = %4.4f\n", tcrit);
    
    for ii = 17 : AddRowExpNodes
        fprintf(fid, '%s\n', Incomplete{ii} );
    end

    % Add shifted Expander Node Coords
    [new_node_9999, NewExpNodes] = P14_UpdateExpanderNodes_V1(original_mesh_nodes, original_mesh_elements, expd_pts_and_normals, offset, theta);
    fprintf(fid,'%d, %6.12f, %6.12f, %6.12f\n',NewExpNodes');

    for ii = AddRowExpNodes+1 : AddRowNode9999
        fprintf( fid, '%s\n', Incomplete{ii} );
    end

    %Add new node 9999 location
    new_node_9999(2:4) = new_node_9999(1:3);
    new_node_9999(1) = 9999;
    fprintf(fid, "%d, %4.6f, %4.6f, %4.6f\n", new_node_9999(1), new_node_9999(2), new_node_9999(3), new_node_9999(4));
   

    for ii = AddRowNode9999+1 : AddRowSpringStiffness
        fprintf( fid, '%s\n', Incomplete{ii} );
    end
    
     
    [Group_1, Group_2, Group_3, Group_4, Group_5, Group_6] = P14_UpdateSx_V1(Sx);

    fprintf(fid, "*Spring, elset=SpringElements4\n");
    fprintf(fid, "\n");
    fprintf(fid, "%4.6f\n", Group_4);
    fprintf(fid, "*Spring, elset=SpringElements3\n");
    fprintf(fid, "\n");
    fprintf(fid, "%4.6f\n", Group_3);
    fprintf(fid, "*Spring, elset=SpringElements6\n");
    fprintf(fid, "\n");
    fprintf(fid, "%4.6f\n", Group_6);
    fprintf(fid, "*Spring, elset=SpringElements2\n");
    fprintf(fid, "\n");
    fprintf(fid, "%4.6f\n", Group_2);
    fprintf(fid, "*Spring, elset=SpringElements5\n");
    fprintf(fid, "\n");
    fprintf(fid, "%4.6f\n", Group_5);
    fprintf(fid, "*Spring, elset=SpringElements1\n");
    fprintf(fid, "\n");
    fprintf(fid, "%4.6f\n", Group_1);


    for ii = AddRowSpringStiffness+1 : AddRowFixed
        fprintf( fid, '%s\n', Incomplete{ii} );
    end
    

    % Add BC
    Fixed_Skin_Nodes=P14_FixSkinNodes_V1(tol, Skin_Nodes, Skin_Elements, NewExpNodes);
    fprintf(fid,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Fixed_Skin_Nodes);
    fprintf(fid, '\n');
    
    for ii = AddRowFixed+1 : length(Incomplete)
        fprintf( fid, '%s\n', Incomplete{ii} );
    end
    fclose( fid );
end

