#!/bin/bash
        # grep  -v '^.............H[^A]' mirrored_1g0y.pdb > mid
        # grep  -v '^............H' mid > mirrored_1g0y_h.pdb
        # grep   '^............HA' mirrored_1g0y.pdb >> mirrored_1g0y_h.pdb
        # sed -i '/HA  GLY/d'  mirrored_1g0y_h.pdb
         tleap -f leap.in
         cat leap.log | grep "Total unperturbed charge:" | awk '{for (i=1;i<=NF;i++);print("1g0y,", $4)}' >> charge.csv  
         mkdir 1g0y
         mv  1g0y_solvated.prmtop 1g0y/cpx_solvated.prmtop
         mv 1g0y_solvated.inpcrd  1g0y/cpx_solvated.inpcrd
         mv 1g0y_solvated.pdb  1g0y/cpx_solvated.pdb
         rm mid leap.in leap.log mirror*.pdb 
         cp -r ctrl_in 1g0y 
         mv md_analysis.sh rmsd.sh rmsf.sh dssp.sh cluster.sh  1g0y
         