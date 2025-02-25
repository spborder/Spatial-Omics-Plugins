"""Codes for generating spot annotations and posting them to an item
"""

import os
import sys

import pandas as pd
import json

from ctk_cli import CLIArgumentParser
import girder_client
import numpy as np

from shapely.geometry import Point
import uuid

#import rpy2.robjects as robjects
import subprocess
from typing_extensions import Union

#from fusion_tools.utils.shapes import load_visium, geojson_to_histomics

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


def load_visium(visium_path:str, include_var_names:list = [], include_obs: list = [], mpp:Union[float,None]=None):
    """Loading 10x Visium Spot annotations from an h5ad file or csv file containing spot center coordinates. Adds any of the variables
    listed in var_names and also the barcodes associated with each spot (if the path is an h5ad file).

    :param visium_path: Path to the h5ad (anndata) formatted Visium data or csv file containing "imagerow" and "imagecol" columns
    :type visium_path: str
    :param include_var_names: List of additional variables to add to the generated annotations (barcode is added by default), defaults to []
    :type include_var_names: list, optional
    :param mpp: If the Microns Per Pixel (MPP) is known for this image then pass it here to save time calculating spot diameter., defaults to None
    :type mpp: Union[float,None], optional
    """

    assert os.path.exists(visium_path)
    """
    if 'h5ad' in visium_path:
        anndata_object = ad.read_h5ad(visium_path)

        if 'spatial' in anndata_object.obsm_keys():

            spot_coords = pd.DataFrame(
                data = anndata_object.obsm['spatial'],
                index = anndata_object.obs_names,
                columns = ['imagecol','imagerow']
            )
        elif all([i in anndata_object.obs_keys() for i in ['imagecol','imagerow']]):
            spot_coords = pd.DataFrame(
                data = {
                    'imagecol': anndata_object.obs['imagecol'],
                    'imagerow': anndata_object.obs['imagerow']
                },
                index = anndata_object.obs_names
            )
    elif 'csv' in visium_path:
        spot_coords = pd.read_csv(visium_path,index_col=0)

    """
    spot_coords = pd.read_csv(visium_path,index_col=0).loc[:,['imagecol','imagerow']]

    # Quick way to calculate how large the radius of each spot should be (minimum distance will be 100um between adjacent spot centroids )
    if mpp is None:
        spot_centers = spot_coords.values
        distance = np.sqrt(
            np.square(
                spot_centers[:,0]-spot_centers[:,0].reshape(-1,1)
            ) + 
            np.square(
                spot_centers[:,1]-spot_centers[:,1].reshape(-1,1)
            )
        )

        min_dist = np.min(distance[distance>0])
        mpp = 1 / (min_dist/100)
    
    # For 55um spot radius
    spot_pixel_radius = int((1/mpp)*55*0.5)

    spot_annotations = {
        'type': 'FeatureCollection',
        'features': [],
        'properties': {
            'name': 'Spots',
            '_id': uuid.uuid4().hex[:24]
        }
    }
    """
    if 'h5ad' in visium_path:
        if len(include_var_names)>0:
            include_vars = [i for i in include_var_names if i in anndata_object.var_names]
        else:
            include_vars = []

        if len(include_obs)>0:
            include_obs = [i for i in include_obs if i in anndata_object.obs_keys()]
        else:
            include_obs = []
    else:
        include_obs = []
        include_vars = []
        
    """
    include_obs = []
    include_vars = []
    anndata_object = None

    barcodes = list(spot_coords.index)
    for idx in range(spot_coords.shape[0]):
        spot = Point(*spot_coords.iloc[idx,:].tolist()).buffer(spot_pixel_radius)

        additional_props = {}
        for i in include_vars:
            additional_props[i] = anndata_object.X[idx,list(anndata_object.var_names).index(i)]
        
        for j in include_obs:
            add_prop = anndata_object.obs[j].iloc[idx]
            if not type(add_prop)==str:
                additional_props[j] = float(add_prop)
            else:
                additional_props[j] = add_prop

        spot_feature = {
            'type': 'Feature',
            'geometry': {
                'type': 'Polygon',
                'coordinates': [list(spot.exterior.coords)]
            },
            'properties': {
                'name': 'Spots',
                '_id': uuid.uuid4().hex[:24],
                '_index': idx,
                'barcode': barcodes[idx]
            } | additional_props
        }

        spot_annotations['features'].append(spot_feature)

    return spot_annotations

def geojson_to_histomics(geojson_anns: Union[list,dict]):

    if type(geojson_anns)==dict:
        geojson_anns = [geojson_anns]
    
    histomics_anns = []
    for g in geojson_anns:
        if 'properties' in g:
            g_name = g['properties']['name']
        else:
            g_name = ''

        histomics_ann = {
            'annotation': {
                'name': g_name,
                'elements': [
                    {
                        'type': 'polyline',
                        'user': f['properties'],
                        'points': [list(i)+[0] if type(i)==tuple else i+[0] for i in f['geometry']['coordinates'][0]]
                    }
                    for f in g['features']
                ]
            }
        }
        histomics_anns.append(histomics_ann)
    
    return histomics_anns


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
    
