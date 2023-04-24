#!/bin/bash
##
## 
##

name=$(basename $0 .sh)
eps=${name}.eps
tmp_dir=_tmp_$name
test -d $tmp_dir || mkdir $tmp_dir

dreg=-R-6/6/-3/3
preg=-R-6.2/6.2/-3.2/3.2
pro=-JX15c/7.5c

inc=.01
dens=50

#grd2stream="../src/grd2stream -r -n 200"
grd2stream="../src/grd2stream"

##
## https://www.wolframalpha.com/input?i=plot+sin%28x*y+%29*cos%28x%29
## create simple topography: z(x,y) = sin(x*y)*cos(x)
## 
gmt grdmath $dreg -I$inc X Y MUL SIN = $tmp_dir/z1.grd
gmt grdmath $dreg -I$inc X COS = $tmp_dir/z2.grd
gmt grdmath $tmp_dir/z1.grd $tmp_dir/z2.grd MUL = $tmp_dir/z.grd

##
## assume flow goes down the hills
##
gmt grdmath $tmp_dir/z.grd DDX -1 MUL = $tmp_dir/u.grd
gmt grdmath $tmp_dir/z.grd DDY -1 MUL = $tmp_dir/v.grd
gmt grdmath $tmp_dir/u.grd $tmp_dir/v.grd R2 SQRT = $tmp_dir/mag.grd

# ncks -O -v z $tmp_dir/u.grd $tmp_dir/u.nc
ncrename -O -v z,u $tmp_dir/u.grd $tmp_dir/u.nc
ncrename -O -v z,v $tmp_dir/v.grd $tmp_dir/v.nc
cdo -O merge $tmp_dir/u.nc $tmp_dir/v.nc  $tmp_dir/velo.nc

## 
## create regular seeds on coarse grid
##
cinc=$(echo "$inc * $dens" | bc)
gmt grdsample $tmp_dir/mag.grd $dreg -I$cinc -G$tmp_dir/mag_c.grd
gmt grdsample $tmp_dir/u.grd $dreg -I$cinc -G$tmp_dir/u_c.grd
gmt grdsample $tmp_dir/v.grd $dreg -I$cinc -G$tmp_dir/v_c.grd
gmt grd2xyz $tmp_dir/mag_c.grd | gmt convert -i0,1 > $tmp_dir/start.d

cat $tmp_dir/start.d | \
  $grd2stream $tmp_dir/velo.nc?u $tmp_dir/velo.nc?v > $tmp_dir/lines_f.d

cat $tmp_dir/start.d | \
  $grd2stream -b $tmp_dir/velo.nc?u $tmp_dir/velo.nc?v > $tmp_dir/lines_b.d


##
## colorscheme
##
cpt=$tmp_dir/topo.cpt
# gmt makecpt -CLakeColor.cpt -T-1/1/.2 > $cpt

## similar color scheme
makecpt -Cocean -T-3/3/.2 > $cpt


########################################################################
### The figure
########################################################################
gmt gmtset PS_MEDIA a4
gmt gmtset PS_PAGE_ORIENTATION portrait
gmt gmtset FONT_TITLE 16p,Helvetica,black

blue=49/98/159  ## uniMS polar blue
red=204/0/0     ## uniMS polar red
green=0/204/0

gmt psbasemap $preg $pro -B0 -K > $eps
gmt grdimage $tmp_dir/z.grd -R -J  -C$cpt -O -K >> $eps

## Topography contours
gmt grdcontour $tmp_dir/z.grd -R -J  -C0.2 -A- -T0.1i/0.05i -O -K >> $eps

## streamlines
test -f $tmp_dir/lines_b.d && \
  gmt psxy -m $tmp_dir/lines_b.d -R -J -W1p,$red -O -K >> $eps

test -f $tmp_dir/lines_f.d && \
  gmt psxy -m $tmp_dir/lines_f.d -R -J -W1p,orange -O -K >> $eps


test -f $tmp_dir/lines.d && \
  gmt psxy -m $tmp_dir/lines.d -R -J -W1p,$green -O -K >> $eps


## -E center on grid node 
## -N do not clip
vec=-Q0.03i/0.1i/0.09in0.25i
vec=-Q0.025i/0.08i/0.06in0.25i
#vec=-Q
gmt grdvector $tmp_dir/u.grd $tmp_dir/v.grd -I0.5 -R -J $vec \
  -G0 -S8i -E -N -O -K  >> $eps

gmt psbasemap -R -J -Ba1f0.2/a1f0.2:."z(x,y) = sin(x*y)*cos(x)":WSne -O >> $eps
echo "%%EOF" >> $eps

pdf=$(basename $eps .eps).pdf
gmt ps2raster -A -Tf $eps

if [ $(uname -s) == "Darwin" ] ; then
  open $pdf
fi

rm -f gmt.conf gmt.history
rm -rf $tmp_dir
