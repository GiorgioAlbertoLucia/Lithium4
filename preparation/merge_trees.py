'''
    Merge multiple trees from a ROOT file into a single tree.
'''

from torchic.utils.terminal_colors import TerminalColors as tc
import pandas as pd
import uproot


if __name__ == '__main__':

    #infile = '/data/galucia/lithium_local/same/LHC24ar_pass2_same_small.root'
    infile = '/data/galucia/lithium_local/same/LHC23_PbPb_pass5_same_small.root'
    table_names = ['O2he3hadtable', 'O2he3hadmult'] #, 'O2he3hadtablemc']
    base = 'DF'

    #outfile = uproot.recreate('/data/galucia/lithium_local/same_merged/LHC24ar_pass2_same_small.root')
    outfile = uproot.recreate('/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass5_same_small.root')

    f = uproot.open(infile)
    keys = list(f.keys())
    _file_folders = [folder for folder in keys if (folder.startswith(base) and '/' not in folder)]
    file_folders_duplicated = [folder.split(';')[0] for folder in _file_folders] # list with potentially duplicated folders
    seen = {}
    for idx, val in enumerate(file_folders_duplicated):
        if val not in seen:
            seen[val] = idx
    file_folders = [_file_folders[idx] for idx in seen.values()]

    for folder in file_folders:
        dfs = []
        for table_name in table_names:
            table_path = f'{infile}:{folder}/{table_name}'
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+table_path+tc.RESET)
            dfs.append(uproot.open(table_path).arrays(library='pd'))

        df = pd.concat(dfs, axis=1)
        df['fIs23'] = True if 'LHC23' in infile else False
        folder_clean = folder.split(';')[0]  # Clean folder name
        outfile[f'{folder_clean}/{table_names[0]}'] = df