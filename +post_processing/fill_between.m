function ff = fill_between(x, y1, y2, rgbcolor, facealpha)

x = vertcat(x, flipud(x));
y = vertcat(y1, flipud(y2));


ff=patch(x,y,rgbcolor,...
  'FaceAlpha', facealpha, 'EdgeColor', 'none');
end