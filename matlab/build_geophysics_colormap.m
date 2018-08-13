close all
j = jet;

%j(1:56,:) ;

figure

p = [1.0 96/255 244/255] ;

%rgb(255, 96, 244)

plot(0,0,'ro','MarkerFaceColor',p,'MarkerEdgeColor',p,...
    'MarkerSize',60);

j(57:72,1) = 1 ;
j(57:72,2) = linspace(0,p(2),16) ;
j(57:72,3) = linspace(0,p(3),16) ;

j(73:88,1) = 1 ;
j(73:88,2) = linspace(p(2),1,16) ;
j(73:88,3) = linspace(p(3),1,16) ;

j = j(1:84,:) ;

peaks; shading interp; view(0,90) ;

% remap to 64 levels
j2=[ interp1(linspace(1,64,84),j(:,1),1:64)', ...
     interp1(linspace(1,64,84),j(:,2),1:64)', ...
     interp1(linspace(1,64,84),j(:,3),1:64)'] ;
colormap(j2); colorbar;




