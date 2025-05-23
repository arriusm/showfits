# showfits.py

sequential fits viewer

with some interactive analysis tools like gauss fit, moffat fit, growing curve analysis and masking

## 
showfits.py V1.4.12 (2025-04-11) (c) USM, Arno Riffeser (arri@usm.lmu.de) based on showfits of Claus Goessl and fitsedit of Johannes
Koppenhoefer

## usage

```
showfits.py [-h] [-figsize FIGSIZE] [-v VERBOSE] [-precaper PRECAPER]
            [-maxaper MAXAPER] [-fitradius FITRADIUS] [-zoom ZOOM]
            [-interpol INTERPOL] [-cmap CMAP] [-full] [-automask]
	    [-autopng] [-autojpg] [-wcslim] [-fixcuts] [-cuts CUTS]
            [imagelist ...]
```

## options
```
  -h, --help            show this help message and exit
  -figsize FIGSIZE      [10.0,8.5] figsize
  -v VERBOSE            [2] verbose
  -precaper PRECAPER    [80] precaper
  -maxaper MAXAPER      [80] maxaper
  -fitradius FITRADIUS  [15] fitradius
  -zoom ZOOM            [0] zoom
  -interpol INTERPOL    [1] interpolation
  -cmap CMAP            [102] cmap
  -full                 [False] full output
  -automask             [False] autosave all masks
  -autopng              [False] autosave all images as png
  -autojpg              [False] autosave all images as jpg
  -wcslim               [False] keep lim according WCS
  -fixcuts              [False] keep cuts fixed from previous image
  -cuts CUTS            [0,0] cuts low and high
```

## usage in plot window
```
          1 - get cursor    at mouse position
          2 - center gauss  at mouse position
          3 - center moffat at mouse position
          4 - growth curve  at mouse position
          5 - psf profile   at mouse position
          6 - star contour  at mouse position
          7 - mask     satellite with two mouse positions
          8 - mask inf satellite with two mouse positions
          0 - mask circle   at mouse position
          b - big mask circle   at mouse position
          m - save mask
          n - save mask
          j - save mask
          @ - UNDO last mask command
          z - decrease lower cut
          Z - increase lower cut
          x - decrease higher cut
          X - increase higher cut
        h r - reset (zoom home)
        < c - zoom back
        > v - zoom forward
          p - pan (move center)
          o - Zoom-to-rect
          s - save image
          f - full frame
          g - Toggle major grids
          G - Toggle minor grids
        k L - x axis log/linear
          l - y axis log/linear
          q - quit
```

## example

```
   > showfits.py  M51*.fits
   nr_images   3
   1   M51_r_stx_220324_103.fits
   2   M51_r_stx_220324_104.fits
   3   M51_r_stx_220324_105.fits
   =============================================================================================================
   nr          1  of  3
   imagename   M51_r_stx_220324_103.fits
   imagesize   3172       4096      
   cuts        1034.0     2192.0    
   cuts        1034.0     1999.0    
   cuts        1034.0     1838.2    
   cuts        1034.0     1704.1    
   cuts        1034.0     1592.4    
   cuts        1034.0     1499.4    
   cuts        1034.0     1421.8    
   cuts        1034.0     1357.2    
   cuts        1034.0     1303.3    
   #                  xc         yc      value
   pixel             564       3262    1027.00
   #                  xc         yc      value
   pixel             956       1979    3281.00
   #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
   gauss          956.15    1978.64    57363.2   1065.910    2503.22       1.60       2.28      76.42       4.50
   #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
   moffat         956.15    1978.64    67269.7   1055.837    2693.50       2.35       3.39      76.51       4.07
   #                  xc         yc    totflux        sky         r0         r1 errtotflux
   grow_curve     956.15    1978.64    63810.5   1060.331       24.0       70.0     1409.5
   grow_curve     956.15    1978.64    63590.8   1060.353       19.0       70.0     1126.3
   grow_curve     956.15    1978.64    63509.3   1060.368       19.0       77.0     1126.3
   ---------------------------------------------------------------------------------------
   grow_curve     956.15    1978.64    63550.2   1060.361       19.0       74.0     1126.3
   ---------------------------------------------------------------------------------------
   ---------------------------------------------------------------------------------------
```