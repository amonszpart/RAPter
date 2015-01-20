#/bin/bash

SVN_ROOT=/export/home/kandinsky/nmellado/svn/globOpt/paper/SIGG2015
SVN_RESULTS=$SVN_ROOT/figures/results

cd bu_eulerCut/
echo "Processing "$PWD

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Euler - Input" --noscatter
python ../normal_distr.py cloudRGBNormal_it45_reProj_noUnass_noPrim.ply ndistrIt45.svg "Euler - Iteration 45"

cp ndistrCloudBinary.svg $SVN_RESULTS/eulercut
cp ndistrIt45.svg $SVN_RESULTS/eulercut

cd ..
cd bu_lHouse3/   
echo "Processing "$PWD

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "House - Input" --noscatter
python ../normal_distr.py cloudRGBNormal_it20_reProj_noUnass_noPrim.ply ndistrIt20.svg "House - Iteration 20"

cp ndistrCloudBinary.svg $SVN_RESULTS/lHouse
cp ndistrIt20.svg $SVN_RESULTS/lHouse

cd ..     
cd bu_empFull/  
echo "Processing "$PWD     

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Empire - Input" --noscatter
python ../normal_distr.py cloudRGBNormal_it10_reProj_noUnass_noPrim.ply ndistrIt10.svg "Empire - Iteration 10"

cp ndistrCloudBinary.svg $SVN_RESULTS/empire
cp ndistrIt10.svg $SVN_RESULTS/empire

cd ..  
cd bu_eulerPM/  
echo "Processing "$PWD       

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "EulerPM - Input" --noscatter
python ../normal_distr.py cloudRGBNormal_it69_reProj_noUnass_noPrim.ply ndistrIt69.svg "Euler - Iteration 69"
python ../normal_distr.py cloudRGBNormal_it69_reProj_noUnass_noPrim.ply ndistrIt69_accum.svg "Euler - Iteration 69 - Acum" --gaussFile2 ../bu_eulerCut/cloudRGBNormal_it45_reProj_noUnass_noPrim.ply 

cp ndistrCloudBinary.svg $SVN_RESULTS/eulerPM
cp ndistrIt69.svg $SVN_RESULTS/eulerPM

cd ..     
cd bu_lansFull/  
echo "Processing "$PWD         

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Lans - Input" --noscatter
python ../normal_distr.py cloudRGBNormal_it32_reProj_noUnass_noPrim.ply ndistrIt32.svg "Lans - Iteration 32"

cp ndistrCloudBinary.svg $SVN_RESULTS/lans
cp ndistrIt32.svg $SVN_RESULTS/lans

cd ..      
cd bu_nolaSigg/   
echo "Processing "$PWD     

python ../normal_distr.py cloud.ply ndistrCloudBinary.svg "Lans - Input" --noscatter
python ../normal_distr.py cloudRGBNormal_it26_reProj_noUnass_noPrim.ply ndistrIt26.svg "Nola - Iteration 26"

cp ndistrCloudBinary.svg $SVN_RESULTS/nola
cp ndistrIt26.svg $SVN_RESULTS/nola

cd ..     
