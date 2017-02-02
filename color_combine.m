function color_combine(r_channel,g_channel)
% r_channel = input('What protein would you like in the red channel? Please type one of the following:\nN-ndc80, ame1, mif2, cep3, ndc10.\n','s');
% b_channel = input('What protein would you like in the blue channel? Please type one of the following:\nN-ndc80, ame1, mif2, cep3, ndc10.\n','s');
rfile = sprintf('heatmap_gray_overlay_%s.png',r_channel);
gfile = sprintf('heatmap_gray_overlay_%s.png',g_channel);
imdata1 = imread(fullfile('\\BioArk.bio.unc.edu\BloomLab\Brandon\Inner Kinetochore Project\Heatmaps\Grayscale Images',rfile));
imdata2 = imread(fullfile('\\BioArk.bio.unc.edu\BloomLab\Brandon\Inner Kinetochore Project\Heatmaps\Grayscale Images',gfile));
imdata3 = imread('\\BioArk.bio.unc.edu\BloomLab\Brandon\Inner Kinetochore Project\Heatmaps\Grayscale Images\heatmap_gray_overlay_template.png');
      RC = imdata1(:,:,1);
      GC = imdata2(:,:,2);
      BC = imdata1(:,:,3) + imdata2(:,:,3);
      im = cat(3,RC,GC,BC);
      figure, imshow(im);
      figure(1);
      imagetitle = sprintf('Magenta Channel: %s | Cyan Channel: %s ',r_channel,g_channel);
      text(280,20,imagetitle,'HorizontalAlignment','center', ...
          'VerticalAlignment','top','FontSize',12,'FontWeight','bold')
end