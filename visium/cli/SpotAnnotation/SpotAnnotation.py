"""Codes for generating spot annotations and posting them to an item
"""

import os
import sys

import pandas as pd
import json

from ctk_cli import CLIArgumentParser
import girder_client

#import rpy2.robjects as robjects
import subprocess

from fusion_tools.utils.shapes import load_visium

INTEGRATION_DATA_KEYS = [
    f'prediction.score.{l}'
    for l in [
        'celltype.l1',
        'celltype.l2',
        'celltype.l3',
        'annotation.l1',
        'annotation.l2',
        'annotation.l3',
        'class',
        'subclass',
        'cluster',
        'cross-species cluster'
    ]
] + ['predsubclassl1','predsubclassl2']


def main(*args):

    gc = girder_client.GirderClient(
        apiUrl=args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments:')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    # Downloading counts file to cwd
    gc.downloadFile(
        args.counts_file,
        path = './'
    )

    file_info = gc.get(f'/file/{args.counts_file}')
    
    # Extracting integration and spot coordinates info
    #extract_spot_info(f'./{file_info["name"]}', INTEGRATION_DATA_KEYS)
    subprocess.call(['Rscript', '../../extract_rds_dataframes.r', *INTEGRATION_DATA_KEYS])

    # Finding all output csv files
    output_csvs = [i for i in os.listdir(os.getcwd()+'/') if 'csv' in i and not i=='spot_coordinates.csv']

    # Creating GeoJSON formatted annotations
    spot_coords = pd.read_csv('spot_coordinates.csv')
    visium_spots = load_visium(spot_coords)

    # Adding properties from other output csv files
    for o in output_csvs:
        property_list = pd.read_csv(o).to_dict('records')
        for s,p in zip(visium_spots['features'],property_list):
            s['properties'] = s['properties'] | p
    
    gc.post(
        f'/annotation/item/{file_info["itemId"]}?token={args.girderToken}',
        data = json.dumps(visium_spots),
        headers = {
            'X-HTTP-Method': 'POST',
            'Content-Type': 'application/json'
        }
    )


if __name__=='__main__':
    main(CLIArgumentParser().parse_args())
    
