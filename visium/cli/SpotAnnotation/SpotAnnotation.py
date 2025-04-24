"""Codes for generating spot annotations and posting them to an item
"""

import os
import sys

import pandas as pd
import json
import geojson

from ctk_cli import CLIArgumentParser
import girder_client
import numpy as np

from shapely.geometry import Point
import uuid

#import rpy2.robjects as robjects
import subprocess
from typing_extensions import Union

from fusion_tools.utils.shapes import load_visium, geojson_to_histomics

# Make sure these don't have spaces
INTEGRATION_DATA_KEYS = [
    f"prediction.score.{l}"
    for l in [
        "celltype.l1",
        "celltype.l2",
        "celltype.l3",
        "annotation.l1",
        "annotation.l2",
        "annotation.l3",
        "class",
        "subclass",
        "cluster",
        "cross-species cluster"
    ]
] + ["predsubclassl1","predsubclassl2"]


def main(args):

    gc = girder_client.GirderClient(
        apiUrl=args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments:')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    file_info = gc.get(f'/file/{args.counts_file}')

    # Downloading counts file to cwd
    gc.downloadFile(
        args.counts_file,
        path = f'{os.getcwd()}/{file_info["name"]}'
    )
    
    # Extracting integration and spot coordinates info
    #extract_spot_info(f'./{file_info["name"]}', INTEGRATION_DATA_KEYS)
    # Sanitizing file name
    file_name_path = f"{os.getcwd()}/{file_info['name']}"
    subprocess.call(['Rscript', '../../extract_rds_dataframes.r', file_name_path,*INTEGRATION_DATA_KEYS])

    if not args.spot_coords is None:
        spot_coords_file_info = gc.get(f'/file/{args.spot_coords}')
        # Downloading spot coordinates file to cwd
        gc.downloadFile(
            args.spot_coords,
            path = f'{os.getcwd()}/spot_coordinates.csv'
        )


    # Finding all output csv files
    output_csvs = [i for i in os.listdir(os.getcwd()+'/') if 'csv' in i and not i=='spot_coordinates.csv']
    print(f'Output CSV files: {output_csvs}')

    # Creating GeoJSON formatted annotations
    spot_coords_path = None
    if 'spot_coordinates.csv' in os.listdir(os.getcwd()+'/'):
        print(f'Spot coordinates found at: {os.getcwd()}/spot_coordinates.csv')
        spot_coords_path = f'{os.getcwd()}/spot_coordinates.csv'
    elif os.path.exists('/spot_coordinates.csv'):
        print(f'Spot coordinates found at: /spot_coordinates.csv')
        spot_coords_path = '/spot_coordinates.csv'
    elif os.path.exists('/cli/spot_coordinates.csv'):
        print(f'Spot coordinates found at: /cli/spot_coordinates.csv')
        spot_coords_path = '/cli/spot_coordinates.csv'

    if not spot_coords_path is None:
        # Loading annotations from spot coordinate path
        visium_spots = load_visium(spot_coords_path)

        # Checking for gene_list_file or gene_selection_method
        if args.use_gene_selection:
            print(f'Using gene selection method: {args.gene_selection_method}, {args.n} selected')
            subprocess.call(['Rscript', '../../gene_selection_csv.r', file_name_path,args.gene_selection_method,str(args.n)])
            output_csvs = [i for i in os.listdir(os.getcwd()+'/') if 'csv' in i and not i=='spot_coordinates.csv']
            print(f'Updated Output CSV files: {output_csvs}')

        if not args.gene_list_file is None:
            try:
                gene_list_file_info = gc.get(f'/file/{args.gene_list_file}')
                print(f'Grabbing specific list of genes from: {gene_list_file_info["name"]}')

                # Downloading counts file to cwd
                gc.downloadFile(
                    args.gene_list_file,
                    path = f'{os.getcwd()}/{gene_list_file_info["name"]}'
                )

                subprocess.call(['Rscript', '../../gene_selection_csv.r', file_name_path,"specific_list",f'{os.getcwd()}/{gene_list_file_info["name"]}'])
                output_csvs = [i for i in os.listdir(os.getcwd()+'/') if 'csv' in i and not i=='spot_coordinates.csv']
                print(f'Updated Output CSV files: {output_csvs}')
            except girder_client.HttpError:
                print('No gene_list_file provided')

        # Adding properties from other output csv files
        for o in output_csvs:
            property_list = pd.read_csv(o).to_dict('records')
            for s,p in zip(visium_spots['features'],property_list):
                s['properties'] = s['properties'] | p
        
        # If a scalefactors_json.json is present
        if not args.scale_factors is None:
            scale_factors_file_info = gc.get(f'/file/{args.scale_factors}')
            gc.downloadFile(
                args.scale_factors,
                path = f'{os.getcwd()}/{scale_factors_file_info["name"]}'
            )

            with open(f'{os.getcwd()}/{scale_factors_file_info["name"]}','r') as f:
                scale_factors = json.load(f)
                f.close()
            
            if 'tissue_hires_scalef' in scale_factors:
                scale_val = scale_factors['tissue_hires_scalef']
            elif 'hires' in scale_factors:
                scale_val = scale_factors['hires']

            visium_spots = geojson.utils.map_geometries(lambda g: geojson.utils.map_tuples(lambda c: (c[0]*scale_val,c[1]*scale_val),g),visium_spots)

        # Converting to histomics format just to add a "name"
        histomics_spots = geojson_to_histomics(visium_spots)

        gc.post(
            f'/annotation/item/{file_info["itemId"]}?token={args.girderToken}',
            data = json.dumps(histomics_spots),
            headers = {
                'X-HTTP-Method': 'POST',
                'Content-Type': 'application/json'
            }
        )


if __name__=='__main__':
    main(CLIArgumentParser().parse_args())
    
