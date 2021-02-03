set term postscript enhanced eps color "Times-Roman" 12
set encoding iso_8859_1
set size ratio 0.8
set lmargin at screen 0.1
#####################################################
set out "Model_SO_PECs.eps"

set title "SO coupled model PECs" font "Times Roman,20"
set xrange [3:16]

set xtics 4,2,16 font "Times Roman,20" 
set xlabel "R ({\305})" font "Times Roman,20"

set yrange [-5500:2000]

set ytics ("-5"-5000,"-4"-4000,"-3"-3000,"-2"-2000,"-1"-1000,"0"0,"1"1000,"2"2000) font "Times Roman,20" 
set ylabel "Energy/hc (1000 cm^{-1})" font "Times Roman,20" offset -3,0

set key samplen 1.5 spacing 2 font ",20" at 19,1000
set grid

set label 1 at 12,1100 "Rb(^2S_{1/2}) + Yb(^3P_2)" font "Times Roman,18"
set label 2 at 12,-500 "Rb(^2S_{1/2}) + Yb(^3P_1)" font "Times Roman,18"
set label 3 at 12,-1300 "Rb(^2S_{1/2}) + Yb(^3P_0)" font "Times Roman,18"

plot "./output_pecs_f90.dat" u 1:2 w l lc rgb '#440154' lw 4 title "1^4{/Symbol P}_{5/2}",\
"./output_pecs_f90.dat" u 1:3 w l lt 2 lw 4 title "2^2{/Symbol P}_{3/2}",\
"./output_pecs_f90.dat" u 1:4 w l lt 3 lw 4 title "1^4{/Symbol P}_{3/2}",\
"./output_pecs_f90.dat" u 1:5 w l lt 4 lw 4 title "1^4{/Symbol S}^+_{3/2}",\
"./output_pecs_f90.dat" u 1:6 w l lc rgb '#fdc328' lw 4 title "2^2{/Symbol P}_{1/2}",\
"./output_pecs_f90.dat" u 1:7 w l lt 6 lw 4 title "1^4{/Symbol P}_{1/2}",\
"./output_pecs_f90.dat" u 1:8 w l lt 7 lw 4 title "1^4{/Symbol P}_{1/2}",\
"./output_pecs_f90.dat" u 1:9 w l lt 8 lw 4 title "3^2{/Symbol S}^+_{1/2}",\
"./output_pecs_f90.dat" u 1:10 w l lt 9 lw 4 title "1^4{/Symbol S}^+_{1/2}"

unset label 1
unset label 2
unset label 3
# plot "./output_pecs.dat" u 1:2 w l lt 8 lw 4 title "3^2{/Symbol S}^+_{1/2}",\
# "./output_pecs.dat" u 1:3 w l lt 6 lw 4 title "1^4{/Symbol P}_{1/2}"

#####################################################
# set out "Model_SO_Diagonal.eps"

# set title "SO matrix diagonal elements" font "Times Roman,18"
# set xrange [3:16]

# set xtics 4,2,16 font "Times Roman,20" 
# set xlabel "R ({\305})" font "Times Roman,20"

# set yrange [-5500:2000]

# set ytics ("-5"-5000,"-4"-4000,"-3"-3000,"-2"-2000,"-1"-1000,"0"0,"1"1000,"2"2000) font "Times Roman,20" 
# set ylabel "Energy/hc (1000 cm^{-1})" font "Times Roman,20" offset -3,0

# set key samplen 1.5 spacing 2 font ",20" at 19,1000
# set grid

# # set label 1 at 12,1100 "Rb(^2S_{1/2}) + Yb(^3P_2)" font "Times Roman,18"
# # set label 2 at 12,-500 "Rb(^2S_{1/2}) + Yb(^3P_1)" font "Times Roman,18"
# # set label 3 at 12,-1300 "Rb(^2S_{1/2}) + Yb(^3P_0)" font "Times Roman,18"

# plot "./diagonal_so-elems.dat" u 1:2 w l lc rgb '#440154' lw 4 title "1^4{/Symbol P}_{5/2}",\
# "./diagonal_so-elems.dat" u 1:3 w l lt 2 lw 4 title "2^2{/Symbol P}_{3/2}",\
# "./diagonal_so-elems.dat" u 1:4 w l lt 3 lw 4 title "1^4{/Symbol P}_{3/2}",\
# "./diagonal_so-elems.dat" u 1:5 w l lt 4 lw 4 title "1^4{/Symbol S}^+_{3/2}",\
# "./diagonal_so-elems.dat" u 1:6 w l lc rgb '#fdc328' lw 4 title "2^2{/Symbol P}_{1/2}",\
# "./diagonal_so-elems.dat" u 1:7 w l lt 6 lw 4 title "1^4{/Symbol P}_{1/2}",\
# "./diagonal_so-elems.dat" u 1:8 w l lt 7 lw 4 title "1^4{/Symbol P}_{1/2}",\
# "./diagonal_so-elems.dat" u 1:9 w l lt 8 lw 4 title "3^2{/Symbol S}^+_{1/2}",\
# "./diagonal_so-elems.dat" u 1:10 w l lt 9 lw 4 title "1^4{/Symbol S}^+_{1/2}"