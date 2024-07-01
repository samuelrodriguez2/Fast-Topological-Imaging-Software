function [Ptes_calcul]=Matrice_imagerie2024(Ptes_calcul,Ptes_expe,Option_normalisation)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function A=Matrice_imagerie(donnees,freq,signal)
% Gradient=A*u
% avec u vecteur des signaux d'excitation du probleme adjoint M(i,omega),
% avec i le ième transducteur
% Gradient est complexe depend des coordonnés du point M. Representer
% (abs(Gradient) correspond à la FTIM
% A dépend DE M, de i, et de omega
% Elle est de taille mxn avec m le nombre de points de l'image et
% n=length(freq)*Nb_transducteurs
% 
% u=t^[M(1,omega(1)) M(1,omega(2)) .... M(1,omega(sf)) , M(2,omega(1)) M(2,omega(2)) .... M(2,omega(sf)) , .... ,  M(N,omega(1)) M(N,omega(2)) .... M(N,omega(sf))  ];
% u(m)=M( E(m/sf)+1  , omega(m-E(m/sf)+1)  )

% A(k,m)=sum( H(i,k,omega(m-E(m/sf)+1)) * S(i,omega(m-E(m/sf)+1)) , i) * H(E(j/sf)+1,k,omega(m-E(m/sf)+1))   


Nb_trans=size(Ptes_expe.fft_signal,2);
Decalage_pas=round(Ptes_expe.Pitch/(Ptes_calcul.X(2)-Ptes_calcul.X(1)));
imag_hors_transducteur_px=round(Nb_trans*0.05)*Decalage_pas;
Larg_imag_px=2*imag_hors_transducteur_px+Nb_trans*Decalage_pas;

% L'abscisse 0 de la fonction de transfert doit être au centre de la source
[~,indice0]=min(abs(Ptes_calcul.X));
indice_bord_image=indice0-imag_hors_transducteur_px;

%% Adding zeros to the transfer functions if image size requires it
ind_min_global=round(indice_bord_image-(Nb_trans-1/2)*Decalage_pas);
ind_max_global=round(indice_bord_image-1/2*Decalage_pas)+Larg_imag_px-1;
if (ind_max_global> size(Ptes_calcul.donnees,2)) || (ind_min_global<1)
    Ptes_calcul.donnees=cat(2,...
        zeros(size(Ptes_calcul.donnees,1), abs(min((ind_min_global-1),0)) ,size(Ptes_calcul.donnees,3)),...
        Ptes_calcul.donnees(:,max(ind_min_global,1):min(ind_max_global,end),:),...
        zeros(size(Ptes_calcul.donnees,1), max(ind_max_global+1-size(Ptes_calcul.donnees,2),0) ,size(Ptes_calcul.donnees,3)) );
    indice0=indice0+abs(min((ind_min_global-1),0));
    indice_bord_image=indice0-imag_hors_transducteur_px;
end
%% Building imaging matrix 
ind_exc=1;
ind_min=round(indice_bord_image-(ind_exc-1/2)*Decalage_pas);
ind_max=ind_min+Larg_imag_px-1;
indices=repmat([ind_min:ind_max]',[1,Nb_trans])-repmat([0:Nb_trans-1]*Decalage_pas,[ind_max-ind_min+1,1]);
indices=reshape(indices,1,[]);
Hi=reshape(Ptes_calcul.donnees(:,indices,:),size(Ptes_calcul.donnees,1),...
    ind_max-ind_min+1,Nb_trans,size(Ptes_calcul.donnees,3));

Hi=reshape(permute(Hi,[3,1,2,4]),Nb_trans,[],size(Ptes_calcul.donnees,3));

Ptes_calcul.X_image=[-imag_hors_transducteur_px:Larg_imag_px-imag_hors_transducteur_px-1]*(Ptes_calcul.X(2)-Ptes_calcul.X(1));
Ptes_calcul.Y_image=Ptes_calcul.Y;
Ptes_calcul.X_trans=([0:Nb_trans-1]+.5)*Ptes_expe.Pitch;


% HS=shiftdim(sum(Hi.*repmat(reshape(Ptes_expe.fft_signal,[Nb_trans 1  length(Ptes_calcul.freq)]),[1 size(Hi,2) 1])  ,1) , 1 );
HS=shiftdim(sum(Hi.*repmat(reshape(Ptes_expe.fft_signal.',[Nb_trans 1  length(Ptes_calcul.freq)]),[1 size(Hi,2) 1])  ,1) , 1 );
if Option_normalisation
    HS=HS./abs(HS).^2;
end
%%juste la transposée de Ptes_expe.fft_signal
Hi=reshape ( shiftdim(Hi,1) , [size(Hi,2), size(Hi,1)*size(Hi,3)]);
Ptes_calcul.A= repmat(HS,1,Nb_trans).*Hi;
