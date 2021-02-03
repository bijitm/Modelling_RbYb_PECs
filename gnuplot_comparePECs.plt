set term postscript enhanced eps color "Times-Roman" 12
set encoding iso_8859_1
set size ratio 0.9
set lmargin at screen 0.1
#####################################################
set out "Compare_PECs.eps"

set title "RbYb: Comparison of ab-initio and model PECs" font "Times Roman,18"
set xrange [3:15]

set xtics 4,2,14 font "Times Roman,20" 
set xlabel "R ({\305})" font "Times Roman,20"

set yrange [-5500:2000]

set ytics ("-5"-5000,"-4"-4000,"-3"-3000,"-2"-2000,"-1"-1000,"0"0,"1"1000,"2"2000) font "Times Roman,20" 
set ylabel "Energy/hc (1000 cm^{-1})" font "Times Roman,20" offset -3,0

set key samplen 0.75 spacing 2 font ",16" at 19.5,1000
set grid

set label 1 at 15.25,1200 "ab initio" font "Times Roman,16"
set label 2 at 17.5,1200 "model" font "Times Roman,16"

plot "./PECs_so.dat" u 1:9 w p lc rgb '#440154' lw 2 title "1^4{/Symbol P}_{5/2}",\
"./PECs_so.dat" u 1:4 w p lt 2 lw 2 title "2^2{/Symbol P}_{3/2}",\
"./PECs_so.dat" u 1:7 w p lt 3 lw 2 title "1^4{/Symbol P}_{3/2}",\
"./PECs_so.dat" u 1:10 w p lt 4 lw 2 title "1^4{/Symbol S}^+_{3/2}",\
"./PECs_so.dat" u 1:2 w p lc rgb '#fdc328' lw 2 title "2^2{/Symbol P}_{1/2}",\
"./PECs_so.dat" u 1:3 w p lt 6 lw 2 title "1^4{/Symbol P}_{1/2}",\
"./PECs_so.dat" u 1:5 w p lt 7 lw 2 title "1^4{/Symbol P}_{1/2}",\
"./PECs_so.dat" u 1:6 w p lt 8 lw 2 title "3^2{/Symbol S}^+_{1/2}",\
"./PECs_so.dat" u 1:8 w p lt 9 lw 2 title "1^4{/Symbol S}^+_{1/2}",\
"./output_pecs_f90.dat" u 1:2 w l lc rgb '#440154' lw 3 title "1^4{/Symbol P}_{5/2}",\
"./output_pecs_f90.dat" u 1:3 w l lt 2 lw 3 title "2^2{/Symbol P}_{3/2}",\
"./output_pecs_f90.dat" u 1:4 w l lt 3 lw 3 title "1^4{/Symbol P}_{3/2}",\
"./output_pecs_f90.dat" u 1:5 w l lt 4 lw 3 title "1^4{/Symbol S}^+_{3/2}",\
"./output_pecs_f90.dat" u 1:6 w l lc rgb '#fdc328' lw 3 title "2^2{/Symbol P}_{1/2}",\
"./output_pecs_f90.dat" u 1:7 w l lt 6 lw 3 title "1^4{/Symbol P}_{1/2}",\
"./output_pecs_f90.dat" u 1:8 w l lt 7 lw 3 title "1^4{/Symbol P}_{1/2}",\
"./output_pecs_f90.dat" u 1:9 w l lt 8 lw 3 title "3^2{/Symbol S}^+_{1/2}",\
"./output_pecs_f90.dat" u 1:10 w l lt 9 lw 3 title "1^4{/Symbol S}^+_{1/2}"

# plot "./PECs_so.dat" u 1:6 w p lt 8 lw 2 title "3^2{/Symbol S}^+_{1/2}",\
# "./PECs_so.dat" u 1:3 w p lt 6 lw 2 title "1^4{/Symbol P}_{1/2}",\
# "./output_pecs.dat" u 1:2 w l lt 8 lw 3 title "3^2{/Symbol S}^+_{1/2}",\
# "./output_pecs.dat" u 1:3 w l lt 6 lw 3 title "1^4{/Symbol P}_{1/2}"

# unset key

#####################################################
# set out "Diagonal_SO_elems.eps"

# set title "Comparison of ab-initio and model diagonal matrix elements" \
#            font "Times Roman,16"
# set xrange [3:15]

# set xtics 4,2,14 font "Times Roman,20" 
# set xlabel "R ({\305})" font "Times Roman,20"

# set yrange [-5500:2000]

# set ytics ("-5"-5000,"-4"-4000,"-3"-3000,"-2"-2000,"-1"-1000,"0"0,"1"1000,"2"2000) font "Times Roman,20" 
# set ylabel "Energy/hc (1000 cm^{-1})" font "Times Roman,20" offset -3,0

# set key samplen 0.75 spacing 2 font ",16" at 19.5,1000
# set grid

# set label 1 at 15.25,1200 "ab initio" font "Times Roman,16"
# set label 2 at 17.5,1200 "model" font "Times Roman,16"

# plot "./PECs_so.dat" u 1:9 w p lc rgb '#440154' lw 2 title "1^4{/Symbol P}_{5/2}",\
# "./PECs_so.dat" u 1:4 w p lt 2 lw 2 title "2^2{/Symbol P}_{3/2}",\
# "./PECs_so.dat" u 1:7 w p lt 3 lw 2 title "1^4{/Symbol P}_{3/2}",\
# "./PECs_so.dat" u 1:10 w p lt 4 lw 2 title "1^4{/Symbol S}^+_{3/2}",\
# "./PECs_so.dat" u 1:2 w p lc rgb '#fdc328' lw 2 title "2^2{/Symbol P}_{1/2}",\
# "./PECs_so.dat" u 1:3 w p lt 6 lw 2 title "1^4{/Symbol P}_{1/2}",\
# "./PECs_so.dat" u 1:5 w p lt 7 lw 2 title "1^4{/Symbol P}_{1/2}",\
# "./PECs_so.dat" u 1:6 w p lt 8 lw 2 title "3^2{/Symbol S}^+_{1/2}",\
# "./PECs_so.dat" u 1:8 w p lt 9 lw 2 title "1^4{/Symbol S}^+_{1/2}",\
# "diagonal_so-elems.dat" u 1:2 w l lc rgb '#440154' lw 3 title "1^4{/Symbol P}_{5/2}",\
# "diagonal_so-elems.dat" u 1:3 w l lt 2 lw 3 title "2^2{/Symbol P}_{3/2}",\
# "diagonal_so-elems.dat" u 1:4 w l lt 3 lw 3 title "1^4{/Symbol P}_{3/2}",\
# "diagonal_so-elems.dat" u 1:5 w l lt 4 lw 3 title "1^4{/Symbol S}^+_{3/2}",\
# "diagonal_so-elems.dat" u 1:6 w l lc rgb '#fdc328' lw 3 title "2^2{/Symbol P}_{1/2}",\
# "diagonal_so-elems.dat" u 1:7 w l lt 6 lw 3 title "1^4{/Symbol P}_{1/2}",\
# "diagonal_so-elems.dat" u 1:8 w l lt 7 lw 3 title "1^4{/Symbol P}_{1/2}",\
# "diagonal_so-elems.dat" u 1:9 w l lt 8 lw 3 title "3^2{/Symbol S}^+_{1/2}",\
# "diagonal_so-elems.dat" u 1:10 w l lt 9 lw 3 title "1^4{/Symbol S}^+_{1/2}"