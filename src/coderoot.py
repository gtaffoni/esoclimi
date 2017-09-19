def coderoot():
 return '/u/gmurante/scratch/Esoclimi/esoMPI-F07/'

def copyfiles(src_dir, template_dir):
 import shutil

 workDir='./'
 shutil.copy(src_dir+"/workarea.py",workDir)
 shutil.copy(src_dir+"/thumblib.py",workDir)
 shutil.copy(src_dir+"/fitslib.py",workDir)
 shutil.copy(src_dir+"/runEBM.py",workDir)
 shutil.copy(src_dir+"/tAtmo.py",workDir)
 shutil.copy(src_dir+"/libraryEBM.py",workDir)
 shutil.copy(src_dir+"/constantsEBM.py",workDir)


