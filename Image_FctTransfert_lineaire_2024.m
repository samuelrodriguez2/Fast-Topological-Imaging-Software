function [X_image,Y_image,ImageTopo,ImageTopoNorm]=Image_FctTransfert_lineaire_2024(Pitch,Ptes_calcul,freq,fft_signal,fft_RF,option_residu)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function [X,Y,ImageTopo,ImageTopoNorm]=Image_FctTransfert_lineaire(Pitch,Ptes_mesure,Ptes_calcul,freq,fft_signal,fft_RF)
%%%%% Entrées
% Pitch : distance inter-élément en m
% Ptes_calcul : structure qui contient la fonction de transfert
% freq : Vecteur qui contient les fréquences de support des variables fft_RF et fft_signal
% fft_signal : Coefficients complexes des transformées de Fourier des signaux d'emission sur les size(fft_signal,1) transducteurs
% fft_RF : Coefficients complexes des transformées de Fourier des signaux mesurés sur les size(fft_RF,1) transducteurs
% option_residu : Si option_residu=1, le residu est retropropage, sinon, c'est uniquement le RF
%%%%% Sorties
% X : Axe des X (parallèle à la barette)
% Y : Axe des Y (normal à la barrette)
% ImageTopo : Matrice contenant les valeurs d'intensité de l'image topologique non normalisée
% ImageTopoNorm : Matrice contenant les valeurs d'intensité de l'image topologique non normalisée par le champ direct

Option_sans_boucle=0;
Nb_exc=size(fft_RF,2);
Decalage_pas=round(Pitch/(Ptes_calcul.X(2)-Ptes_calcul.X(1)));

if abs(Decalage_pas-Pitch/(Ptes_calcul.X(2)-Ptes_calcul.X(1)))>1e-4
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('Attention, le Pitch n''est pas completement compatible avec le pas d''espace de la fonction de transfert ! \n');
    fprintf('Rapport entre le pitch et la pas d''espace : %2.6f\n',Pitch/(Ptes_calcul.X(2)-Ptes_calcul.X(1)));
    fprintf('Risque d''image distordue\n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
end

[val,indice0]=min(abs(Ptes_calcul.X));
Decalage_pas=round(Pitch/(Ptes_calcul.X(2)-Ptes_calcul.X(1)));
imag_hors_transducteur_px=round(Nb_exc*0.05)*Decalage_pas;
Larg_imag_px=2*imag_hors_transducteur_px+Nb_exc*Decalage_pas;
indice_bord_image=indice0-imag_hors_transducteur_px;

%%%%%%%%%%%% Calcul des champs %%%%%%%%%%%%%%%%%%%%%
spectre_obj1=ones(1,1,size(fft_signal,1),'single'); % Spectre pour un point prÃ©cis
spectre_obj2=ones(1,1,size(fft_RF,1),'single'); % Spectre pour un point prÃ©cis
champ_direct= zeros(size(Ptes_calcul.donnees,1),Larg_imag_px,size(Ptes_calcul.donnees,3),'single');
champ_adjoint=zeros(size(Ptes_calcul.donnees,1),Larg_imag_px,size(Ptes_calcul.donnees,3),'single');
%% Adding zeros to the transfer functions if image size requires it
ind_min_global=round(indice_bord_image-(Nb_exc-1/2)*Decalage_pas);
ind_max_global=round(indice_bord_image-1/2*Decalage_pas)+Larg_imag_px-1;
if (ind_max_global> size(Ptes_calcul.donnees,2)) || (ind_min_global<1)
    abs(min((ind_min_global-1),0)) 
    max(ind_max_global-size(Ptes_calcul.donnees,2),0)
    Ptes_calcul.donnees=cat(2,...
        zeros(size(Ptes_calcul.donnees,1), abs(min((ind_min_global-1),0)) ,size(Ptes_calcul.donnees,3)),...
        Ptes_calcul.donnees(:,max(ind_min_global,1):min(ind_max_global,end),:),...
        zeros(size(Ptes_calcul.donnees,1), max(ind_max_global+1-size(Ptes_calcul.donnees,2),0) ,size(Ptes_calcul.donnees,3)) );
    indice0=indice0+abs(min((ind_min_global-1),0));
    indice_bord_image=indice0-imag_hors_transducteur_px;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if option_residu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [val,indice0]=min(abs(Ptes_calcul.X));
    if mod(indice0,2)==0
        donnees_gamma=(Ptes_calcul.donnees(1,2:Decalage_pas:end,:));
        X_gamma=Ptes_calcul.X(2:Decalage_pas:end);
    else
        donnees_gamma=(Ptes_calcul.donnees(1,1:Decalage_pas:end,:));
        X_gamma=Ptes_calcul.X(1:Decalage_pas:end);
    end
    [val,indice0]=min(abs(X_gamma));
    donnees_gamma(1,indice0-3:indice0+3,round(end/2));
    u_gamma= zeros(size(donnees_gamma,1),Nb_exc,size(donnees_gamma,3),'single');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calcul de la solution du problème direct au niveau des transducteurs
    for ind_exc=1:Nb_exc
        spectre_obj1(1,1,:)=(single(fft_signal(:,ind_exc)));
        ind_min=indice0-ind_exc+1;
        ind_max=ind_min+Nb_exc-1;
        
        if (ind_max> size(donnees_gamma,2)) || (ind_min<1)
            donnees_recalees=cat(2,...
                zeros(size(donnees_gamma,1), abs(min((ind_min-1),0)) ,size(donnees_gamma,3)),...
                donnees_gamma(:,max(ind_min,1):min(ind_max,end),:),...
                zeros(size(donnees_gamma,1), max(ind_max-size(donnees_gamma,2),0) ,size(donnees_gamma,3)) );
        else
            donnees_recalees=donnees_gamma(:,ind_min:ind_max,:);
        end
        u_gamma=u_gamma  +repmat(spectre_obj1,size(donnees_recalees,1),size(donnees_recalees,2)).*donnees_recalees;
    end
    u_gamma=(shiftdim(u_gamma,1));    
    fft_RF=fft_RF-u_gamma.';
end

if Option_sans_boucle
    ind_min=round(indice_bord_image-1/2*Decalage_pas);
    ind_max=ind_min+Larg_imag_px-1;
    indices=repmat([ind_min:ind_max]',[1,Nb_exc])-repmat([0:Nb_exc-1]*Decalage_pas,[ind_max-ind_min+1,1]);
    indices=reshape(indices,1,[]);
    H=reshape(Ptes_calcul.donnees(:,indices,:),size(Ptes_calcul.donnees,1),ind_max-ind_min+1,Nb_exc,[]);
    
    champ_direct=  squeeze(sum(H.*repmat(shiftdim(fft_signal.',-2),[size(Ptes_calcul.donnees,1),ind_max-ind_min+1]),3));
    champ_adjoint= squeeze(sum(H.*repmat(shiftdim(fft_RF',-2),[size(Ptes_calcul.donnees,1),ind_max-ind_min+1]),3));
else
    for ind_exc=1:Nb_exc
        % Attention : image inversÃ©e par rapport Ã  d'habitude
        spectre_obj1(1,1,:)=single(fft_signal(:,ind_exc));
        spectre_obj2(1,1,:)=single(conj(fft_RF(:,ind_exc)));
        ind_min=round(indice_bord_image-(ind_exc-1/2)*Decalage_pas);
        ind_max=ind_min+Larg_imag_px-1;
        
        donnees_recalees=Ptes_calcul.donnees(:,ind_min:ind_max,:);
        
        champ_direct =champ_direct +repmat(spectre_obj1,size(donnees_recalees,1),size(donnees_recalees,2)).*donnees_recalees;
        champ_adjoint=champ_adjoint+repmat(spectre_obj2,size(donnees_recalees,1),size(donnees_recalees,2)).*donnees_recalees;
    end
end

X_image=[-imag_hors_transducteur_px:Larg_imag_px-imag_hors_transducteur_px-1]*(Ptes_calcul.X(2)-Ptes_calcul.X(1));
Y_image=Ptes_calcul.Y;

ImageTopo=squeeze((sum(champ_direct.*champ_adjoint,3)));
ImageTopo=ImageTopo/(max(max(abs(ImageTopo))));
ImageTopoNorm=squeeze((sum(champ_direct.*champ_adjoint./abs(champ_direct).^2,3)));
ImageTopoNorm=ImageTopoNorm/(max(max(abs(ImageTopoNorm))));

