import os
import shutil
from warnings import *
from Bio import PDB, BiopythonWarning
from pyrosetta import *


simplefilter('ignore', BiopythonWarning)
pdb_io = PDB.PDBIO()
pdb = PDB.PDBParser()
pdb_name=[]
for num_id, pdbfile in enumerate([file for file in os.listdir() if file.endswith('.pdb')]):
    num_id = num_id + 1  
    decoy_dir = os.path.basename(pdbfile).split('.pdb')[0]
    pdb_name.append(decoy_dir)
    pdb_name.sort()

def make_leap(name):
    with open('leap.in', 'w') as f:
        line = ''' source leaprc.protein.ff14SB
        source leaprc.water.tip3p
        pro = loadpdb %s.pdb
        #saveamberparm protein %s.prmtop %s.inpcrd
        solvatebox  pro TIP3PBOX 12
        charge pro
        addionsrand pro Na+ 0
        addionsrand pro Cl- 0
        savepdb pro %s_solvated.pdb
        saveamberparm pro %s_solvated.prmtop %s_solvated.inpcrd
        quit
        ''' % (name,name,name,name,name,name)
        f.write(line)

def make_md(name):
    with open('md_analysis.sh', 'w') as f:
        line = '''#!/bin/bash %s
pmemd.cuda -O -i ctrl_in/min1.in -p cpx_solvated.prmtop -o min1.out -r min1.rst -ref cpx_solvated.inpcrd -c cpx_solvated.inpcrd
pmemd.cuda -O -i ctrl_in/min2.in -p cpx_solvated.prmtop -o min2.out -r min2.rst -c min1.rst
pmemd.cuda -O -i ctrl_in/heat.in -o heat.out -p cpx_solvated.prmtop -c min2.rst -r heat.rst -x heat.mdcrd -ref min2.rst
pmemd.cuda -O -i ctrl_in/density.in -o density.out -p cpx_solvated.prmtop -c heat.rst -r density.rst -x density.mdcrd -ref heat.rst
pmemd.cuda -O -i ctrl_in/equil.in -o equil.out -p cpx_solvated.prmtop -c density.rst -r  equil.rst -x equil.mdcrd -ref density.rst
pmemd.cuda -O -i ctrl_in/prod.in -o prod.out -p cpx_solvated.prmtop -c equil.rst -r prod.rst -x prod.mdcrd -ref equil.rst
#analysis
bash rmsd.sh
bash rmsf.sh
bash cluster.sh
bash dssp.sh
sed '12 aset terminal png truecolor' -i dssp.gnu
sed '13 aset output "dssp.png"' -i dssp.gnu
gpuplot dssp.gnu &
        ''' %(name)
        f.write(line)

def rmsd_sh(write_path):
    txt = '''\
cpptraj <<EOF
parm cpx_solvated.prmtop
trajin prod.mdcrd 1 last 1
reference cpx_solvated.pdb
rms reference out rmsd.txt :1-20@CA,C,O,N&!@H=
EOF
'''
    with open(os.path.join(write_path, 'rmsd.sh'), 'w') as f:
        f.write(txt)

def rmsf_sh(write_path):
    txt = '''\
cpptraj <<EOF
parm cpx_solvated.prmtop
trajin prod.mdcrd
rms first
average crdset MyAvg
run
rms ref MyAvg
atomicfluct out rmsf.dat @C,CA,N byres
EOF
'''
    with open(os.path.join(write_path, 'rmsf.sh'), 'w') as f:
        f.write(txt)

def dssp_sh(write_path):
    txt = '''\
cpptraj <<EOF
parm cpx_solvated.prmtop
trajin prod.mdcrd
secstruct :1-20 out dssp.gnu sumout dssp.agr
EOF
'''
    with open(os.path.join(write_path, 'dssp.sh'), 'w') as f:
        f.write(txt)

def cluster_sh(write_path):
    txt = '''\
cpptraj <<EOF
parm cpx_solvated.prmtop
trajin prod.mdcrd
strip :WAT
strip :Na+
cluster C0 \\
        hieragglo clusters 5 averagelinkage \\
        rms :1-20@C,N,O,CA,CB&!@H=  \\
        out cnumvtime.dat \\
        sil Sil \\
        summary summary.dat \\
        info info.dat \\
        cpopvtime cpopvtime.agr normframe \\
        repout rep repfmt pdb \\
        singlerepout singlerep.nc singlerepfmt netcdf \\
        avgout Avg avgfmt restart
EOF
'''
    with open(os.path.join(write_path, 'cluster.sh'), 'w') as f:
        f.write(txt)

def make_sh(name):
    with open('run.sh', 'w') as f:
         line = '''#!/bin/bash
        # grep  -v '^.............H[^A]' mirrored_%s.pdb > mid
        # grep  -v '^............H' mid > mirrored_%s_h.pdb
        # grep   '^............HA' mirrored_%s.pdb >> mirrored_%s_h.pdb
        # sed -i '/HA  GLY/d'  mirrored_%s_h.pdb
         tleap -f leap.in
         cat leap.log | grep "Total unperturbed charge:" | awk '{for (i=1;i<=NF;i++);print("%s,", $4)}' >> charge.csv  
         mkdir %s
         mv  %s_solvated.prmtop %s/cpx_solvated.prmtop
         mv %s_solvated.inpcrd  %s/cpx_solvated.inpcrd
         mv %s_solvated.pdb  %s/cpx_solvated.pdb
         rm mid leap.in leap.log mirror*.pdb 
         cp -r ctrl_in %s 
         mv md_analysis.sh rmsd.sh rmsf.sh dssp.sh cluster.sh  %s
         ''' %(name,name,name,name,name,name,name,name,name,name,name,name,name,name,name)
         f.write(line)

for i in pdb_name :
#   m = MirrorPDBMover(i+str('.pdb'))
#    m.apply()
    make_leap(i)
    make_md(i)
    rmsd_sh('.')
    rmsf_sh('.')
    dssp_sh('.')
    cluster_sh('.')
    make_sh(i)
    os.system('sh run.sh')

#    os.mkdir(i)
#    os.chdir(str(i))
#    shutil.move(str(i)+'_solvated.prmtop',str(i)/'cpx_solvated.prmtop')
#    shutil.move(str(i)+'_solvated.inpcrd',str(i)/'cpx_solvated.inpcrd')
#    os.chdir('../')
