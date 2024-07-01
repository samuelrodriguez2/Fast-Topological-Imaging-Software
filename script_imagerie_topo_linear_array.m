%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2024 Samuel Rodriguez
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choosing the computation method
Option_imaging=1;
% 1: First generation FTIM     G=UxV         [S. Rodriguez et al., Ultrasonics 52, 1010–1018, 2012, 10.1016/j.ultras.2012.08.002]
% 2: Fourier-based simulation  G=omega^2xUxV [S. Rodriguez et al., Ultrasonics 54, 1880–1890, 2014, 10.1016/j.ultras.2013.10.001]
% 3: matricial FTIM            G=UxV         [S. Rodriguez et al., Ultrasonics 68, 8–16, 2016, 10.1016/j.ultras.2016.02.002]
% NB: Option 3 is more RAM-demanding
%% Loading experimental data
load('Aluminum_piece_with_7holes.mat');
% Single plane wave insonification

Screen_size=get(0,'ScreenSize');
figure_size=round([1 Screen_size(4)*.15 [0.95 0.75].*Screen_size(3:4)]);

figure('Position',figure_size);
h1=subplot(132);
imagesc([1,size(Ptes_expe.signalIn,2)],[0 , (size(Ptes_expe.signalIn,1)-1)/Ptes_expe.fe*1e6],Ptes_expe.signalIn,[-1 1])
xlabel('Transducer');ylabel('Time [µs]');title('Excitation signals')
set(gca,'FontSize',18)
colormap jet;colorbar
hold on

h2=subplot(133);
imagesc([1,size(RFdata,2)],[0 , (size(RFdata,1)-1)/Ptes_expe.fe*1e6],RFdata)
xlabel('Transducer');ylabel('Time [µs]');title('RF')
set(gca,'FontSize',18)
colormap jet;colorbar
hold on

subplot(131)
[X,map] = imread('experimental_setup.jpg');
imshow(X,map)
title('Experimental setup');set(gca,'FontSize',18)


drawnow

%% Defining the array and the reference medium used for imaging
% Ptes_expe.tStart corresponds to the spatial offset.
% It corrects the impedance adaptation layer of the array
% tStart and space-offset relation is given by:
% tStart=2*offset/C;

% Offset truncation
ind_tstart=round(Ptes_expe.tStart*Ptes_expe.fe)+1;
RF=RFdata(ind_tstart:end,:);
Ptes_expe.N_trans=size(RF,2);

% Image
Ptes_expe.Profondeur=60e-3; % Image depth
Pitch_divider=2;            % Needs to be an integer
Ptes_expe.resolution=Ptes_expe.Pitch/Pitch_divider;
Ptes_expe.tRange=round(2.2*Ptes_expe.Profondeur / Ptes_expe.C * Ptes_expe.fe) /Ptes_expe.fe;
Ptes_expe.t = [0:1/Ptes_expe.fe:(Ptes_expe.tRange)-1/Ptes_expe.fe];

%% Preprocessing the experimental data before imaging
% Setting the saturated RF signal beginning to 0
T_erase=3.2e-6; % [µs]
[~,ind_T_erase]=min(abs(Ptes_expe.t-T_erase));
RF(1:ind_T_erase,:)=0;
plot([1,size(RF,2)],[1 1]*T_erase*1e6,'LineWidth',2,'Color',[1 0 0],'LineStyle','--','Parent',h2);

% Harmonization between the measurement duration and the spatial depth of the image
if length(Ptes_expe.t)>size(RF,1)
    warning('The required depth is greater than the corresponding time duration of the RF\n RF will be zero-padded\n');
    RF=[RF;zeros(length(Ptes_expe.t)-size(RF,1),size(RF,2))];
else
    RF=RF(1:length(Ptes_expe.t),:);
    plot([1,size(RF,2)],[1 1]*Ptes_expe.t(end)*1e6,'LineWidth',2,'Color',[1 0 0],'LineStyle','--','Parent',h2);
end

if length(Ptes_expe.t)>size(Ptes_expe.signalIn,1)
    Ptes_expe.signalIn=[Ptes_expe.signalIn;zeros(length(Ptes_expe.t)-size(Ptes_expe.signalIn,1),size(Ptes_expe.signalIn,2))];
else
    Ptes_expe.signalIn=Ptes_expe.signalIn(1:length(Ptes_expe.t),:);
    plot([1,size(Ptes_expe.signalIn,2)],[1 1]*Ptes_expe.t(end)*1e6,'LineWidth',2,'Color',[1 0 0],'LineStyle','--','Parent',h1);
end
drawnow

% Computation of frequency-domain signals
Ptes_expe.freq=[0:size(Ptes_expe.signalIn,1)-1]/size(Ptes_expe.signalIn,1)*Ptes_expe.fe;
Ptes_expe.fft_signal=fft(Ptes_expe.signalIn);
fft_RF = fft(RF);

% Frequency range of interest
BF=Ptes_expe.f0*[0.5 1.5]; 
[~,ind_fmin]=min(abs(Ptes_expe.freq-BF(1)));
[~,ind_fmax]=min(abs(Ptes_expe.freq-BF(2)));

Ptes_calcul.freq=Ptes_expe.freq(ind_fmin:ind_fmax); 
Ptes_expe.fft_signal=Ptes_expe.fft_signal(ind_fmin:ind_fmax,:);
fft_RF = fft_RF(ind_fmin:ind_fmax,:);


%% Computing the images
if Option_imaging~=2
    fprintf('Computing the radiation transfer function of one element : ')
    tic
    Ptes_calcul=FT_transducteur(Ptes_expe.C,Ptes_expe.Largeur_transducteur,round(Ptes_expe.N_trans*2.2)*Ptes_expe.Largeur_transducteur,Ptes_expe.Profondeur,Ptes_expe.resolution,Ptes_calcul.freq);
    h=toc; fprintf('%f s\n',h)
end

switch Option_imaging
    case 1
        fprintf('Computing directly the image (FTIM): ')
        tic
        [Ptes_calcul.X_image,Ptes_calcul.Y_image,ImageTopo,ImageTopoNorm]=Image_FctTransfert_lineaire_2024(Ptes_expe.Pitch,Ptes_calcul,Ptes_expe.freq,Ptes_expe.fft_signal,fft_RF,0);
        h=toc; fprintf('%f s\n',h)
    case 2
        fprintf('Computing using Fourier propagation \n')
        tic;
        [Ptes_calcul.X_image,Ptes_calcul.Y_image,ImageTopo,ImageTopoNorm]=...
            Image_Fourier_isotrope(Ptes_expe.Pitch,Pitch_divider,[0:Ptes_expe.resolution:Ptes_expe.Profondeur],Ptes_calcul.freq.',Ptes_expe.fft_signal,fft_RF,Ptes_expe.C,1.1);
        h=toc; fprintf('%f s\n',h)
    case 3
        fprintf('Computing imaging matrix (used in SelF-EASE) : ')
        Option_normalisation=0;
        tic;        
        [Ptes_calcul]=Matrice_imagerie2024(Ptes_calcul,Ptes_expe,Option_normalisation);
        h=toc; fprintf('%f s\n',h)
        
        % Construction de l'image
        fprintf('Computing the image using the matrix: ')
        tic;
        u = reshape(conj(single(fft_RF)),size(fft_RF,1)*size(fft_RF,2),1);
        ImageTopo_ligne=(Ptes_calcul.A*u); %%G=AX, cela forme l'image du milieu
        ImageTopo=reshape(ImageTopo_ligne,length(Ptes_calcul.Y_image),[]);%size(Ptes_calcul.A,1)/length(Ptes_calcul.Y_image));
        ImageTopo=ImageTopo/(max(max(abs(ImageTopo))));
        h=toc; fprintf('%f s\n',h)
end
%% Plots
Echelle=[0,1]*max(max(abs(ImageTopo)));
Bornes_X=[Ptes_calcul.X_image(1) Ptes_calcul.X_image(end)];
Bornes_Y=[Ptes_calcul.Y_image(1) Ptes_calcul.Y_image(end)];

figure_size2=figure_size+round([0.2*Screen_size(3) 0 0 0]);
figure('Position',figure_size2);
if Option_imaging~=3
    subplot(121)
    imagesc(Ptes_calcul.X_image*1e3,1e3*Ptes_calcul.Y_image,abs(ImageTopo),Echelle);
    axis image;
    set(gca,'XLim',Bornes_X*1e3,'YLim',Bornes_Y*1e3,'FontSize',18);
    axis ij;
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    title('without normalization')
    colormap jet
    colorbar
    
    subplot(122)
    imagesc(Ptes_calcul.X_image*1e3,1e3*Ptes_calcul.Y_image,abs(ImageTopoNorm));
    axis image;
    set(gca,'XLim',Bornes_X*1e3,'YLim',Bornes_Y*1e3,'FontSize',18);
    axis ij;
    xlabel('x axis [mm]')
    % ylabel('y axis [mm]')
    title('with normalization')
    colormap jet
    colorbar
else
    imagesc(Ptes_calcul.X_image*1e3,1e3*Ptes_calcul.Y_image,abs(ImageTopo),Echelle);
    axis image;
    set(gca,'XLim',Bornes_X*1e3,'YLim',Bornes_Y*1e3,'FontSize',18);
    axis ij;
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    colormap jet
    colorbar
end
