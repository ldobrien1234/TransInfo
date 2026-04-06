function make_plots(datafile)
%This function will make several plots of TI curves given the data file's
%path relative to the working directory
%Args:
%   datafile (string): file path relative to working directory
%Returns:
%   None

    load(datafile);

    figure;
    imshow(imdata);
    hold on
    scatter((pxc+bxc),(pyc-byc),'c*');
    hold on
    scatter(round(N/2),round(0.94*M),300,'r.');
    hold off
    
    figure;
    surf(transRange,scaleRange,TI);
    xlabel('Translating');
    ylabel('Scaling');
    
    figure;
    surf(transRange,scaleRange,TI);
    xlabel('Translating');
    ylabel('Scaling');
    xlim([0,200])
    ylim([0.8,1])
    
    
    % TI=imgaussfilt(TI,2);
    [minimaMask,~] = islocalmax2(-TI,'MinProminence',0.1);
    [maximaMask,~] = islocalmax2(TI,'MinProminence',0.1);
    
    %Extract coordinates
    [X,Y]=meshgrid(transRange,scaleRange);
    
    % Plot
    figure;
    contour(transRange, scaleRange, TI, 50); 
    hold on;
    plot3(X(minimaMask),Y(minimaMask),TI(minimaMask), 'r.');
    hold on
    plot3(X(maximaMask),Y(maximaMask),TI(maximaMask),'b.');
    title('Approximate Symmetries');
    xlim([0,800])
    ylim([0.5,1])
    xlabel('Translating');
    ylabel('Scaling');
    legend('Contours', 'Symmetries');
    hold off;

end