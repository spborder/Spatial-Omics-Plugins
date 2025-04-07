"""Processing uploaded Xenium data
"""

import os

from fusion_tools.utils.shapes import (
    load_polygon_csv,
    align_object_props,
    geojson_to_histomics
)
from ctk_cli import CLIArgumentParser

import girder_client
import pandas as pd
import json


def main(args):

    gc = girder_client.GirderClient(
        apiUrl = args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments:')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    # Input arguments:
    # coords_file: Coordinates for either cell boundaries or cell centroids
    # cell_info: Cell-level data where one column is "cell-id"
    # image_item: Item to which generated annotations will be attached

    coords_file_info = gc.get(f'/file/{args.coords_file}')
    # Downloading coords and cell-info files
    gc.downloadFile(
        args.coords_file,
        path = f'{os.getcwd()}/{coords_file_info["name"]}'
    )

    if not args.cell_info is None:
        cell_file_info = gc.get(f'/file/{args.cell_info}')
        gc.downloadFile(
            args.cell_info,
            path = f'{os.getcwd()}/{cell_file_info["name"]}'
        )
        cell_file = pd.read_csv(f'{os.getcwd()}/{cell_file_info["name"]}')


    # Getting just the column names in csv file
    coords_file_columns = pd.read_csv(f'{os.getcwd()}/{coords_file_info["name"]}',nrows=0).columns.tolist()
    print(f'Coordinates file columns: {coords_file_columns}')

    if all([i in coords_file_columns for i in ['x_centroid','y_centroid']]):
        # This file contains just centroid information
        shape_cols = ['x_centroid','y_centroid']
        shape_options = {
            'radius': 'cell_area' if 'cell_area' in coords_file_columns else 'nucleus_area'
        }


    elif all([i in coords_file_columns for i in ['x_vertex','y_vertex']]):
        # This file contains cell boundary coordinates
        shape_cols = ['x_vertex','y_vertex']
        shape_options = {}

    cell_annotations = load_polygon_csv(
        csv_path = f'{os.getcwd()}/{coords_file_info["name"]}',
        name = 'Xenium Cells',
        shape_cols = shape_cols
    )

    if not args.cell_info is None:
        cell_annotations = align_object_props(
            cell_annotations,
            cell_file,
            prop_cols = [i for i in cell_file if not i=='cell_id'],
            alignment_type = 'cell_id'
        )

    cell_annotations_histomics = geojson_to_histomics(cell_annotations)

    gc.post(
        f'/annotation/item/{args.image_item}',
        data = json.dumps(cell_annotations_histomics)
    )



if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

