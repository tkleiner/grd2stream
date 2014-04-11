
# This code is only demostration code it does not work as an script

lonMin=-98
lonMax=-71
latMin=18
latMax=33

# Create basemap instance
map = Basemap(projection='merc',
              llcrnrlat=latmin,
              urcrnrlat=latmax,
              llcrnrlon=lonmin,urcrnrlon=lonmax,resolution="f"
              ,fix_aspect=True)


# lon,lat are the coordinate variable

# Change coordinates to Basemap Projection
x, y =  map(lon, lat)

# plot the background map
lines = 5;
latStep = abs((latmax-latmin)/lines);
lonStep = abs((lonmax-lonmin)/lines);
map.drawparallels(np.arange(latmin+latStep,latmin+(lines)*latStep,latStep),
                  labels=[0,1,0,0],color='black',fontsize=7,zorder=20,
                  linewidth=0.8,fmt=gf.lat2str,rotation=270);
map.drawmeridians(np.arange(lonmin+lonStep,lonmin+(lines)*lonStep,lonStep),
                  labels=[0,0,1,0],color='black',fontsize=7,zorder=20,
                  linewidth=0.8,fmt=gf.lon2str)

map.drawcoastlines(linewidth=0.5,color='#000000');
map.drawmapboundary();
map.fillcontinents();



# Data is from the u,v components of the current and c is the module c=sqrt(v*v+u*u)


plt.streamplot(x,y,u,v,density=(7,7), color="#000000", linewidth=3*c/c.max());


Q = plt.pcolormesh(x,y,c,shading='gouraud',alpha=1, cmap=mpl.cm.Blues);

# the colorbar parameters are set to 0,2 min,max and the shrink to 0.5
Q.set_clim(vmin=barMin, vmax=barMax);
cb = colorbar(shrink=barShrink);
cb.ax.set_title('Current speed \n knots',fontsize=8);

