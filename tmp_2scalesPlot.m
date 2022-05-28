
 label_font = 10;
 left = 0.1; bottom = 0.1; delta = 0.06 ; %delta = 0.03
 width = 0.2 ; height = 0.5  ;
 vshift = height + delta ; hshift = width + delta ; 

figure(fig_number) ; clf 
   gray = [ 0.7 0.7 0.7 ] ;
   blue = [ 0.0 1.0 1.0 ] ;
   % hleft = axes('Position', [left bottom width height ]) ;
   %hcent = axes('Position', [left+hshift  bottom width height ]) ;
   
   hplt1 = plot( localizationList, mseFnLoc(:,1)/(sigma2_l + sigma2_s), ...
                 'k' , 'LineWidth', 2) ; hold on
   hplt2 = plot( localizationList, mseFnLoc(:,2)/(sigma2_l + sigma2_s), ...
                 'b' , 'LineWidth', 1.25) ; hold on
   hplt3 = plot( localizationList, mseFnLoc(:,3)/(sigma2_l + sigma2_s), ...
                 'r' , 'LineWidth', 1.25) ; hold on
             
   hplt1 = plot( localizationList, mseTransfFnLoc(:,1)/(sigma2_l + sigma2_s), ...
                 'k:' , 'LineWidth', 2) ; hold on
   hplt2 = plot( localizationList, mseTransfFnLoc(:,2)/(sigma2_l + sigma2_s), ...
                 'b:' , 'LineWidth', 1.25) ; hold on
   hplt3 = plot( localizationList, mseTransfFnLoc(:,3)/(sigma2_l + sigma2_s), ...
                 'r:' , 'LineWidth', 1.25) ; hold on
             
   axis_save = axis ; axis_save(3) = 0; axis( axis_save ) ;
   title(paramString,'FontSize',2)
   xlabel( 'R (grid points)', 'FontSize', label_font ) ; 
   ylabel( 'MSE / (prior MSE)', 'FontSize', label_font ) ; 
     %set( gca,'XTickLabel','||||||||||')
     % set( gca,'YTick', [0:2:14] )
   set( gca, 'FontSize',label_font, 'TickDir', 'out') ;
   
   figure(fig_number+1) ; clf % "production figure" ----------------------
   gray = [ 0.7 0.7 0.7 ] ;
   blue = [ 0.0 1.0 1.0 ] ;
   % hleft = axes('Position', [left bottom width height ]) ;
   %hcent = axes('Position', [left+hshift  bottom width height ]) ;
   
   h = axes('Position', [.2 .2 .7 .7] )
   hplt1 = plot( localizationList, mseFnLoc(:,1)/(sigma2_l + sigma2_s), ...
                 'k' , 'LineWidth', 2) ; hold on
   hplt2 = plot( localizationList, mseFnLoc(:,2)/(sigma2_l + sigma2_s), ...
                 'b' , 'LineWidth', 2) ; hold on
   hplt3 = plot( localizationList, mseFnLoc(:,3)/(sigma2_l + sigma2_s), ...
                 'r' , 'LineWidth', 2) ; hold on
             
   hplt1 = plot( localizationList, mseTransfFnLoc(:,1)/(sigma2_l + sigma2_s), ...
                 'k:' , 'LineWidth', 2) ; hold on
   hplt2 = plot( localizationList, mseTransfFnLoc(:,2)/(sigma2_l + sigma2_s), ...
                 'b:' , 'LineWidth', 2) ; hold on
   hplt3 = plot( localizationList, mseTransfFnLoc(:,3)/(sigma2_l + sigma2_s), ...
                 'r:' , 'LineWidth', 2) ; hold on
             
   axis( [0 40 0 0.5] ); axis square
   label_font = 12 ;
   xlabel( 'localization cutoff (grid points)', 'FontSize', label_font ) ; 
   ylabel( 'MSE / (prior MSE)', 'FontSize', label_font ) ; 
     %set( gca,'XTickLabel','||||||||||')
    set( gca,'YTick', [0:0.1:0.5] ) ; set( gca,'XTick', [0:10:40] )
   set( gca, 'FontSize',label_font, 'TickDir', 'out') ;
exportfig( gcf, ['figExample_test2.eps'], 'Width', 6.5, 'Color', 'rgb' )
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hcent = axes('Position', [left+hshift  bottom width height ]) ;
   hplt1 = plot( localizationList, mse_l(4,:)./mse_l(2,:), ...
                 'Color', blue , 'LineWidth', 2) ; hold on
                 % 'b' , 'LineWidth', 2) ; hold on
   hplt2 = plot( localizationList, mse_l(3,:)./mse_l(2,:), ...
                 'Color', blue , 'LineWidth', 0.5) ; 
                % 'b' , 'LineWidth', 0.5) ; 
   hplt3 = plot( localizationList, mse_l(1,:)./mse_l(2,:), ...
                'b--', 'Color', blue , 'LineWidth', 0.25) ; 
                 % 'b:'                ) ; 
   xlabel( 'R (grid points)', 'FontSize', label_font ) ; 
   set( gca, 'FontSize',label_font, 'TickDir', 'out') ;
   
   hrght = axes('Position', [left+2*hshift  bottom width height ]) ;
   hplt1 = plot( localizationList, mse(4,:)./mse(2,:), ...
                 'k' , 'LineWidth', 2) ; hold on
   hplt2 = plot( localizationList, mse(3,:)./mse(2,:), ...
                 'k' , 'LineWidth', 0.5) ; 
   hplt3 = plot( localizationList, mse(1,:)./mse(2,:), ...
                 'k--'               ) ; 
   xlabel( 'R (grid points)', 'FontSize', label_font ) ; 
   set( gca, 'FontSize',label_font, 'TickDir', 'out') ;
   axis_save = axis, axis_save(3) = 0; axis( axis_save ) ;

exportfig( gcf, ['fig_' expt_name '.eps'], 'Width', 6.5, 'Color', 'cmyk' )
