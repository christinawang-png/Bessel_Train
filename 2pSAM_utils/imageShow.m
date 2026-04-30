function imageShow(imageSction)
 %%%

  %newimageSction = log(imageSction + eps(0));
  b = min(min(imageSction));
  fprintf(num2str(b));
  fprintf('\n');
  a = max(max(imageSction));
  fprintf(num2str(a));
  fprintf('\n');
  %b = min(min(newimageSction));
  %a = max(max(newimageSction));


  k = 1.0/(a-b);
  scaleImageSction = k*(imageSction-b);

  colormap(hot); 
  set(gca,'ColorScale','log'); 
  imagesc(log(imageSction + eps(0))); 
  %clim([-8 log(a)]); 
  clim([-7 0]); 
  axis equal; 
  axis off; 
  drawnow;
end