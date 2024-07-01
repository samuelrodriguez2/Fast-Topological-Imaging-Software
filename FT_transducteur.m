function Ptes_calcul=FT_transducteur(C,Largeur_transducteur,Largeur,Profondeur,Pas,freq)
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
% function Ptes_calcul=FT_transducteur(C,Largeur_transducteur,Largeur,Profondeur,Pas,freq)
% C : wave velocity
% Largeur_transducteur : transducer width
% Largeur : Width of the computed radiation pattern 
% Profondeur : Depth of the computed radiation pattern 
% Pas : spatial step
% freq : frequency vector
Lambda_min=C/max(freq);
Ptes_calcul.freq=freq;
Pas_integration=1/(ceil(Largeur_transducteur/2/(Lambda_min/10))/(Largeur_transducteur/2));
Nb_intervalles_int=Largeur_transducteur/Pas_integration;
%% Discrétisation
Xmin=-round(Largeur/2/Pas)*Pas;
Xmax=Xmin+Largeur;

Ptes_calcul.X=[Xmin:Pas:Xmax];
Ptes_calcul.Y=[Pas:Pas:Profondeur];
omega=2*pi*shiftdim(Ptes_calcul.freq,-1);

[X,Y] = meshgrid(Ptes_calcul.X,Ptes_calcul.Y);
x0=[0.5:Nb_intervalles_int/2-.5]*Pas_integration;

%% Computation of the radiation pattern
Ptes_calcul.donnees=zeros([size(X),length(Ptes_calcul.freq)],'single');
for ind=1:length(x0)    
    R=single(sqrt((X-x0(ind)).^2+Y.^2));
    Beta=single(repmat(omega/C,size(R)));
%     Ptes_calcul.donnees=Ptes_calcul.donnees+1i*besselh(0,Beta.*repmat(R,size(omega)));
    Ptes_calcul.donnees=Ptes_calcul.donnees+single(1i*conj(besselh(0,Beta.*repmat(R,size(omega)))));
end

Ptes_calcul.donnees=(Ptes_calcul.donnees+Ptes_calcul.donnees(:,end:-1:1,:))*Pas_integration;


