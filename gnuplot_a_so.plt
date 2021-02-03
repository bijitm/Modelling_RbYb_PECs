set term postscript enhanced eps color "Times-Roman" 12
set encoding iso_8859_1
# set size 1,1
set size ratio 0.75
#####################################################
set out "a_so_func.eps"
# set view 60,75
# set hidden3d
set title "R-dependent a_{so,Yb}" font "Times Roman,20"
set xrange [3:16]

set xtics 4,2,16 font "Times Roman,20" 
set xlabel "R ({\305})" font "Times Roman,20"

# set yrange [450:850]

# set ytics 450,50,850 font "Times Roman,20" 
set ytics font "Times Roman,20" 
set ylabel "a_{so,Yb}(R) (cm^{-1})" font "Times Roman,20" offset -3,0

set key spacing 2 font ",20" #bottom
set grid

# set label 1 at 10,620 "a_{so,Yb}exp(-{/Symbol a}(R-R_0)^4)" font "Times Roman,20"
# set label 2 at 10,580 "a_{so,Yb} = 806.09 cm^{-1}" font "Times Roman,20"
# set label 3 at 10,540 "{/Symbol a} = 0.00002 {\305}^{-4}" font "Times Roman,20"
# set label 4 at 10,500 "R_0 = 16 {\305}" font "Times Roman,20"
# set label 5 at 10,660 "Guess function:" font "Times Roman,20"


# plot "./soc_function.dat" u 1:2 w l lw 4 title "a_{so,Yb}(R) guess",\
# "./output_least_sq_fit.dat" u 1:6 w l lw 4 title "a_{so,Yb}(R) fit"
# plot "./soc_function.dat" u 1:2 w l lw 4 title "a_{so,Yb}(R)"

plot "./soc_function_f90.dat" u 1:2 w l lw 4 title "a_1(R)",\
"./soc_function_f90.dat" u 1:3 w l lw 4 title "a_2(R)",\
"./soc_function_f90.dat" u 1:4 w l lw 4 title "a_3(R)",\
"./soc_function_f90.dat" u 1:5 w l lw 4 title "a_4(R)"

