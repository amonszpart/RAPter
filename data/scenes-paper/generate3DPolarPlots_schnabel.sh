#/bin/bash

SVN_ROOT=/export/home/kandinsky/nmellado/svn/globOpt/paper/SIGG2015
SVN_COMPARAISONS=$SVN_ROOT/figures/comparisons/schnabel


cd bu_lHouse3/   
echo "Processing "$PWD

python ../normal_distr.py schnabel_minsup10_pts.ply ndistr_schnabel_minsup10.svg "lHouse - Schnabel (minsup=10)"

cp ndistr_schnabel_minsup10.svg $SVN_COMPARAISONS/lHouse

cd ..     
cd bu_empFull/  
echo "Processing "$PWD     

python ../normal_distr.py schnabel_minsup20000_pts.ply ndistr_schnabel_minsup20000.svg "Lans - Schnabel (minsup=20000)"

cp ndistr_schnabel_minsup20000.svg $SVN_COMPARAISONS/emp

cd ..  
cd bu_eulerPM/  
echo "Processing "$PWD       

python ../normal_distr.py schnabel_minsup50_pts.ply ndistr_schnabel_minsup50.svg "eulerPM - Schnabel (minsup=50)"

cp ndistr_schnabel_minsup50.svg $SVN_COMPARAISONS/eulerPM

cd ..     
cd bu_lansFull/  
echo "Processing "$PWD         

python ../normal_distr.py schnabel_minsup1000_pts.ply ndistr_schnabel_minsup1000.svg "Lans - Schnabel (minsup=1000)"

cp ndistr_schnabel_minsup1000.svg $SVN_COMPARAISONS/lans

cd ..      
cd bu_nolaSigg/   
echo "Processing "$PWD     

python ../normal_distr.py schnabel_minsup500_pts.ply ndistr_schnabel_minsup500.svg "Lans - Schnabel (minsup=500)"

cp ndistr_schnabel_minsup500.svg $SVN_COMPARAISONS/nola

cd ..     
