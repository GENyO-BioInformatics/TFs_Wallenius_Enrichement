import os, glob
# venv/bin/python3 analyse_random_lists.py --tfs data/TFsLists/WikiPathways/size3/169_TFsAnnoted.txt

python = '/home/adriangarciamoreno/Desktop/TFs_Wallenius_Enrichement/venv/bin/python3'
TFsListsFiles = glob.glob('data/TFsLists/*/*/*_TFsAnnoted.txt'); len(TFsListsFiles)

# Select the annotation to use
# total = len(TFsListsFiles)
cmds = []; i = 1
TFsListFile = TFsListsFiles[0]
for TFsListFile in TFsListsFiles:
    outfile = TFsListFile.replace('TFsLists','EnrResults').replace('.txt','.tsv')
    if not os.path.exists(outfile):
      #print('Doing {} of {} for {}'.format(i,total, annotation))
      cmd = '{} analyse_random_lists.py --tfs {} # {}'.format(python,TFsListFile, i)
      cmds.append(cmd); i += 1

with open('null_simulations.sh','w') as cmdsFile:
    cmdsFile.write('\n'.join(cmds))


