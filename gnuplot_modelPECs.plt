set term postscript enhanced eps color "Times-Roman" 12
set encoding iso_8859_1
# set size 1,1
set size ratio 0.75
#####################################################
set out "Model_PECs.eps"
# set view 60,75
# set hidden3d
set title "RbYb spin-free model PECs after fitting" font "Times Roman,20"
set xrange [3:16]

set xtics 4,2,16 font "Times Roman,20" 
set xlabel "R ({\305})" font "Times Roman,20"

set yrange [-5500:2000]

set ytics ("-5"-5000,"-4"-4000,"-3"-3000,"-2"-2000,"-1"-1000,"0"0,"1"1000,"2"2000) font "Times Roman,20" 
set ylabel "Energy/hc (1000 cm^{-1})" font "Times Roman,20" offset -3,0

set key spacing 2 font ",20" bottom invert
set grid

set label 1 at 12.75,-350 "Rb(^2S) + Yb(^3P)" font "Times Roman,20"

plot "./model_pecs.dat" u 1:2 w l lw 4 title "3^2{/Symbol S}^+ (V_1)",\
"./model_pecs.dat" u 1:3 w l lw 4 title "2^2{/Symbol P} (V_2)",\
"./model_pecs.dat" u 1:4 w l lw 4 title "1^4{/Symbol P} (V_3)",\
"./model_pecs.dat" u 1:5 w l lw 4 title "1^4{/Symbol S}^+ (V_4)"


#####################################################
# set out "Model_PECs_python.eps"
# # set view 60,75
# # set hidden3d
# set title "RbYb spin-free guess PECs before fitting" font "Times Roman,20"
# set xrange [3:16]

# set xtics 4,2,16 font "Times Roman,20" 
# set xlabel "R ({\305})" font "Times Roman,20"

# set yrange [-5500:2000]

# set ytics ("-5"-5000,"-4"-4000,"-3"-3000,"-2"-2000,"-1"-1000,"0"0,"1"1000,"2"2000) font "Times Roman,20" 
# set ylabel "Energy/hc (1000 cm^{-1})" font "Times Roman,20" offset -3,0

# set key spacing 2 font ",20" bottom invert
# set grid

# set label 1 at 12.75,-350 "Rb(^2S) + Yb(^3P)" font "Times Roman,20"

# # plot "./model_pecs.dat" u 1:2 w l lw 4 title "3^2{/Symbol S}^+ (V_1)",\
# # "./model_pecs.dat" u 1:4 w l lw 4 title "1^4{/Symbol P} (V_3)"

# plot "./model_pecs.dat" u 1:2 w l lw 4 title "3^2{/Symbol S}^+ (V_1)",\
# "./model_pecs.dat" u 1:3 w l lw 4 title "2^2{/Symbol P} (V_2)",\
# "./model_pecs.dat" u 1:4 w l lw 4 title "1^4{/Symbol P} (V_3)",\
# "./model_pecs.dat" u 1:5 w l lw 4 title "1^4{/Symbol S}^+ (V_4)"

# # plot "./model_pecs_python.dat" u 1:2 w l lw 4 title "3^2{/Symbol S}^+ (V_1)",\
# # "./model_pecs_python.dat" u 1:3 w l lw 4 title "2^2{/Symbol P} (V_2)",\
# # "./model_pecs_python.dat" u 1:4 w l lw 4 title "1^4{/Symbol P} (V_3)",\
# # "./model_pecs_python.dat" u 1:5 w l lw 4 title "1^4{/Symbol S}^+ (V_4)"