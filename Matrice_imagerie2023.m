function Ptes_calcul=Matrice_imagerie2023(Ptes_calcul,Ptes_expe)
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
% donnees_lin=reshape(Ptes_calcul.donnees,size(Ptes_calcul.donnees,1)*size(Ptes_calcul.donnees,2),size(Ptes_calcul.donnees,3));
Decalage_pas=round(Ptes_expe.Pitch/(Ptes_calcul.X(2)-Ptes_calcul.X(1)));
imag_hors_transducteur_px=round(Nb_trans*0.05)*Decalage_pas;
Larg_imag_px=2*imag_hors_transducteur_px+Nb_trans*Decalage_pas;
%Larg_imag_px=(Nb_trans*Decalage_pas);

%%%Hi matrice de -' 128 x 147 x 81
Hi=zeros(Nb_trans,size(Ptes_calcul.donnees,1)*Larg_imag_px,size(Ptes_calcul.donnees,3),'single');
HS=zeros(1,size(Hi,2),size(Hi,3),'single');

% L'abscisse 0 de la fonction de transfert doit être au centre de la source
[~,indice0]=min(abs(Ptes_calcul.X));
indice_born_image=indice0-imag_hors_transducteur_px;
% donnees_recalees=zeros(size(Ptes_calcul.donnees,1),Larg_imag_px,size(Ptes_calcul.donnees,3),'single');

for ind_exc=1:Nb_trans
    % ind_min et ind_max sont les bornes que l'on va prendre dans la
    % fonction de transfert pour le transducteur ind_exc
    ind_min=round(indice_born_image-(ind_exc-1/2)*Decalage_pas);
    ind_max=ind_min+Larg_imag_px-1;

    if (ind_max> size(Ptes_calcul.donnees,2)) || (ind_min<1)
        fprintf('Matrice_imagerie : La fonction de transfert %d est complétée par des zeros pour atteindre la largeur désirée\n',ind_exc)
        donnees_recalees=cat(2,...
            zeros(size(Ptes_calcul.donnees,1), abs(min((ind_min-1),0)) ,size(Ptes_calcul.donnees,3)),...
            Ptes_calcul.donnees(:,max(ind_min,1):min(ind_max,end),:),...
            zeros(size(Ptes_calcul.donnees,1), max(ind_max-size(Ptes_calcul.donnees,2),0) ,size(Ptes_calcul.donnees,3)) );
    else
        donnees_recalees=Ptes_calcul.donnees(:,ind_min:ind_max,:);
    end
    Hi(ind_exc,:,:)=shiftdim(reshape(donnees_recalees,size(donnees_recalees,1)*size(donnees_recalees,2),size(donnees_recalees,3)),-1);
    HS=HS+Hi(ind_exc,:,:).*repmat(shiftdim(Ptes_expe.fft_signal(:,ind_exc).',-1),[1 size(Hi,2) 1]);
end
HS=shiftdim(HS,1);    

clear('donnees_recalees');

Ptes_calcul.X_image=[-imag_hors_transducteur_px:Larg_imag_px-imag_hors_transducteur_px-1]*(Ptes_calcul.X(2)-Ptes_calcul.X(1));
Ptes_calcul.Y_image=Ptes_calcul.Y;
Ptes_calcul.X_trans=([0:Nb_trans-1]+.5)*Ptes_expe.Pitch;


Hi2=reshape ( shiftdim(Hi,1) , [size(Hi,2), size(Hi,1)*size(Hi,3)]);
clear('Hi');

Ptes_calcul.A= repmat(HS,1,Nb_trans).*Hi2;

