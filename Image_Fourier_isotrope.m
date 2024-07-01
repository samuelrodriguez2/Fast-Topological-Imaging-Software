function [CoorReel,zz,ImageTopo,ImageTopoNorm]=Image_Fourier_isotrope(Pitch,DivX,zz,freq,fft_signal,fft_RF,vPhase,Nb_max)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(Pitch,DivX,BF,zz,Ptes_mesure,RF,vPhase,Nb_max)

if isempty(DivX)
    Pitch0=Pitch;
else
    Pitch0=Pitch;
    Pitch=Pitch/DivX;
    fft_RF2=zeros(size(fft_RF,1),size(fft_RF,2)*DivX);
    fft_signal2=zeros(size(fft_signal,1),size(fft_signal,2)*DivX);
    for ind=1:DivX
        fft_RF2(:,ind:DivX:end-DivX+ind)=fft_RF;
        fft_signal2(:,ind:DivX:end-DivX+ind)=fft_signal;
    end
    fft_RF=fft_RF2;
    fft_signal=fft_signal2;
end
Nb_exc=size(fft_RF,2);    

CoorReel=single(Pitch*[0:Nb_exc-1]);
if size(fft_RF,2)==1;CoorReel=CoorReel+Pitch/2;end;

%%%%%%%%%%%%%%%%%%%%%
z=zeros(1,1,length(zz));
z(1,1,:)=zz;

%%%%%%%%% Vecteur et vitesses de phase
if min(size(vPhase))==2 % si vecteur frequence vPhase
    vPhase=single(interp1(vPhase(:,1),vPhase(:,2),freq,'pchip'));
else % si polynome ou cste
    vPhase=single(polyval(vPhase,freq));
end
% fig=figure;plot(freq*1e-6,vPhase,'-');
% set(gca,'FontSize',16)
% xlabel('Frequency [MHz]')
% ylabel('Phase Velocity [m/s]')
% drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Agrandissement de la barette (pour éviter des éffets de repliements selon x)
% Cet agrandissement dépend de la directivité des transducteur cad de ka
% Mal Calculé !!!, il faudrait prendre en compte le k_min de l'ensemble de
% la barette
ka_min=min(2*pi*freq./vPhase)*Pitch0;
Ltotale=CoorReel(end)-CoorReel(1)+Pitch;
duree_totale=1/diff(freq(1:2));
if ka_min<5
    elargissement=max(vPhase)*duree_totale/2;
else

    elargissement=max(vPhase)*duree_totale/4;
end
elargissement=min(elargissement,Pitch*0.5*(Nb_max/(length(zz)*length(freq))-length(CoorReel)));

elargissement=zz(end);
fprintf('Taille de l''elargissement latéral anti-repliement : %3.3f \n',elargissement);

vec1=[-Pitch:-Pitch:-elargissement];vec1=vec1(end:-1:1);
vec2=[CoorReel(end)+Pitch:Pitch:elargissement+CoorReel(end)];
Coor= [ vec1 , CoorReel , vec2 ];

Ltotale=Coor(end)-Coor(1)+Pitch;
kx=single(2*pi*[0:round(1/Pitch*Ltotale-1)]/Ltotale);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=length(freq);
fft_ligne_adjoint=fft([zeros(l1,length(vec1),'single')   conj(fft_RF)   zeros(l1,length(vec2),'single')],[],2);
% 

ligne_direct=[zeros(l1,length(vec1),'single') fft_signal zeros(l1,length(vec2),'single')];
fft_ligne_direct= fft(ligne_direct,[],2);


%%%%%%%% Matrice depropagation delon z %%%%%%%%%%%%
l2=size(fft_ligne_direct,2);
l3=length(z);
propag=zeros(l1,l2,l3,'single');
kx_corr=[kx(1:round(length(kx)/2)) kx(round(length(kx)/2)+1:end)-kx(end)-kx(2)];


kz=sqrt(repmat(((2*pi*freq)./(vPhase)).^2,[1,l2,l3]) - repmat((kx_corr).^2,[l1,1,l3]) );
propag=exp(-1i*repmat(z,[l1,l2,1]).* real(kz))-(real(kz)==0); % Terme de propagation nul si kz imaginaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Calcul des champs %%%%%%%%%%%%%%%%%%%%%
% avec passage de l'expace (kx,z) à l'espace (x,z)
fft_champ_direct=repmat(fft_ligne_direct,[1,1,l3]).*propag;
fft_champ_direct=fft_champ_direct./repmat(((2*pi*freq)./(vPhase)),[1,l2,l3]).*kz;

champ_direct=ifft(fft_champ_direct,[],2);
% champ_direct=champ_direct./repmat(((2*pi*freq)./(vPhase)),[1,l2,l3]).*kz;

clear fft_champ_direct;
ind_x=[length(vec1)+1 length(vec1)+Nb_exc];
champ_direct =champ_direct(:,ind_x(1):ind_x(2),:);

fft_champ_adjoint=repmat(fft_ligne_adjoint,[1,1,l3]).*propag;
%%%%%%%%%%% division par cn
fft_champ_adjoint=fft_champ_adjoint.*(repmat(2*pi*freq./vPhase,[1,l2,l3])./kz);
clear propag kz;
champ_adjoint=ifft(fft_champ_adjoint,[],2);
clear fft_champ_adjoint;
champ_adjoint=champ_adjoint(:,ind_x(1):ind_x(2),:);

mat_freq= repmat(freq,[1,size(champ_direct,2),size(champ_direct,3)]).^(2);
ImageTopo=    squeeze((sum((2*pi*mat_freq).*champ_direct.*champ_adjoint,1)));
ImageTopoNorm=squeeze((sum((2*pi*mat_freq).*champ_direct.*champ_adjoint./abs(champ_direct).^2,1)));


ImageTopo=ImageTopo.'/(max(max(abs(ImageTopo))));
ImageTopoNorm=ImageTopoNorm.'/(max(max(abs(ImageTopoNorm))));

