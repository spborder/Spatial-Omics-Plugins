"""Implementing cell composition deconvolution
"""
import os
import sys

from ctk_cli import CLIArgumentParser
import girder_client

#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri
import subprocess


ORGAN_REF_KEY = {
    "Azimuth Adipose Reference": "adiposeref",
    "Azimuth Bone Marrow Reference": "bonemarrowref",
    "Azimuth Fetus Reference": "fetusref",
    "Azimuth Heart Reference": "heartref",
    "Azimuth Human Cortex Reference": "humancortexref",
    "Azimuth Kidney Reference": "kidneyref",
    "KPMP Atlas Kidney": "kidneykpmp",
    "Azimuth Lung Reference": "lungref",
    "Azimuth Pancreas Reference": "pancreasref",
    "Azimuth Mouse Cortex Reference": "mousecortexref",
    "Azimuth PBMC Reference": "pbmcref",
    "Azimuth Tonsil Reference": "tonsilref",
    "ST Deconvolve": "st_deconvolve"
}


INTEGRATION_DATA_KEYS = {
    'adiposeref': ['celltype.l1','celltype.l2'],
    'bonemarrowref': ['celltype.l1','celltype.l2'],
    'fetusref': ['annotation.l2','annotation.l1'],
    'heartref': ['celltype.l1','celltype.l2'],
    'humancortexref': ['class','subclass','cluster','cross-species cluster'],
    'kidneyref': ['annotation.l1','annotation.l2','annotation.l3'],
    'kidneykpmp': ['subclass.l2','subclass.l1'],
    'lungref': ['annotation.l1','annotation.l2'],
    'pancreasref': ['celltype.l1','celltype.l2'],
    'mousecortexref': ['class','subclass','cluster','cross-species cluster'],
    'pbmcref': ['celltype.l1','celltype.l2','celltype.l3'],
    'tonsilref': ['celltype.l1','celltype.l2']
}

# pancreas might be in there also ['annotation.l1']
# liver might be in there also ['celltype.l1','celltype.l2']
# tonsil has a v2, ['celltype.l1','celltype.l2']
# mouse pansci = ['Main_cell_type']
# lung has a v1 and v2, v2 = ['ann_level_1','ann_level_2','ann_level_3','ann_level_4','ann_level_5','ann_finest_level']


def main(args):

    gc = girder_client.GirderClient(
        apiUrl = args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments:')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    if not args.organ == 'Not Listed':
        # print contents of current working directory, see if files were copied over
        print('Contents of working directory')
        print(os.listdir(os.getcwd()+'/'))
        file_info = gc.get(f'/file/{args.counts_file}')

        # Downloading counts file to cwd
        gc.downloadFile(
            args.counts_file,
            path = f'{os.getcwd()}/{file_info["name"]}'
        )
        print('Updated contents of directory')
        print(os.listdir(os.getcwd()+'/'))

        print(f'Running cell deconvolution for: {args.organ}')
        #integrator = robjects.globalenv['integrate_spatial']
        #integrator(
        #    file_info['name'],
        #    ORGAN_REF_KEY[args.organ]
        #)

        subprocess.call(['Rscript', '.././.cell_deconvolution.r', '"'+file_info['name']+'"','"'+ORGAN_REF_KEY[args.organ]+'"'])

        print(os.listdir(os.getcwd()+'/'))

        print(f'Uploading file to {file_info["itemId"]}')
        # Posting integration results to item
        file_ext = file_info['name'].split('.')[-1]
        uploaded_file = gc.uploadFileToItem(
            itemId = file_info['itemId'],
            filepath = f'./{file_info["name"].replace("."+file_ext,"_integrated.rds")}'
        )
        
if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

