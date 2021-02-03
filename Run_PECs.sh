set -e
f95 PECs.f90 -o PECs -L/home/bijit/Softwares/lapack-3.9.0/ -llapack -lrefblas && ./PECs
gnuplot "./gnuplot_modelPECs.plt"
gnuplot "./gnuplot_abinitioPECs.plt"
gnuplot "./gnuplot_modelSOPECs.plt"
gnuplot "./gnuplot_comparePECs.plt"
gnuplot "./gnuplot_a_so.plt"
cp *.eps "../Documents_results/"
cd "../Documents_results"
latex "Spin_orbit_coupling_RbYb.tex" >log
dvipdf "Spin_orbit_coupling_RbYb" 2> log
cd -