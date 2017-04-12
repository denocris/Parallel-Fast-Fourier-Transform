set term png
set output "diffusivity.png"
minT=0
maxT=0.9
set cbrange [minT:maxT]
plot 'diffusivity0002.dat' matrix with image

set term gif animate
set output "animate.gif"
frames = 4 
minT=0
maxT=0.02
set cbrange [minT:maxT]
do for [i=1:frames] {
  if (i<10){
     plot 'concentration000'.i.'.dat' matrix  with image
     }
  if (i>=10  && i<100){
     plot 'concentration00'.i.'.dat' matrix  with image
     }
  if (i>= 100 && i<1000){
     plot 'concentration0'.i.'.dat' matrix  with image
     }
  if (i>= 1000 && i<10000){
     plot 'concentration'.i.'.dat' matrix  with image
     }
}

