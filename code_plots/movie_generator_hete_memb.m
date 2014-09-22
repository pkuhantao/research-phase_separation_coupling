%% movie generator

function movie_generator_hete_memb
close all;
clear all;
global folderNAME frameNUM frameINDEX frameSTEP num_frames_per_second

caseNUM = 4; %% total number of cases
caseINDEX = 0; %% staring index of case
caseSTEP = 1; %% case step

num_frames_per_second = 10;  % output frequency

%% loop for different cases
for p=1:1:caseNUM
    index = (p-1)*caseSTEP + caseINDEX;

folderNAME = strcat('/Users/than/Desktop/research_data/simulation_cps/Yin_Yang_grids/Sep_20_2014/cps/hete_mismatch/phi_s_zero_psi_m_nonzero/case',...
            num2str(index));

frameNUM = 78;  % total number of frames
frameINDEX = 0; %% starting index of frame
frameSTEP = 200000; %% frame step

%% movie for inner solvent and the membrane
generateMovie('/CPS_pics/insol_psi_minor_phase.mp4', '/CPS_pics/inSolv_YinYang_psi_minor_phase_');
generateMovie('/CPS_pics/insol_psi_solid.mp4', '/CPS_pics/inSolv_YinYang_psi_solid_');
generateMovie('/CPS_pics/memb_psi.mp4', '/CPS_pics/memb_YinYang_psi_');

end

end

%% movie generationg function
function generateMovie(movieNAME, imageNAME) 
global folderNAME frameNUM frameINDEX frameSTEP num_frames_per_second

mvOBJ = VideoWriter( strcat(folderNAME, movieNAME), 'MPEG-4');
mvOBJ.FrameRate = num_frames_per_second;
open(mvOBJ);
for p=1:1:frameNUM
    index = (p-1)*frameSTEP + frameINDEX;
    fileNAME = strcat(folderNAME, imageNAME, num2str(index), '.png');
    I = imread(fileNAME, 'png');
    frame = im2frame(I);    
    writeVideo(mvOBJ, frame)
end
%  Once all the frames have been generated and added, we must close the file.
close(mvOBJ);
end
