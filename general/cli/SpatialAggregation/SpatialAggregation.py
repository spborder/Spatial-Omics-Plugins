"""Plugin for performing spatial aggregation from one set of annotations (aggregator) to another/others
"""
import os
import sys
import numpy as np
import json

import girder_client
from ctk_cli import CLIArgumentParser

from fusion_tools.handler.dsa_handler import DSAHandler
from fusion_tools.utils.shapes import spatially_aggregate, export_annotations


def main(args):

    sys.stdout.flush()

    # Initialize girder client
    gc = girder_client.GirderClient(
        apiUrl = args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments: ')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    file_info = gc.get(f'/file/{args.input_image}')
    image_id = file_info['itemId']
    image_name = file_info['name']

    dsa_handler = DSAHandler(
        girderApiUrl=args.girderApiUlr
    )
    annotations = dsa_handler.get_annotations(
        item = image_id
    )

    ann_names = [i['properties']['name'] for i in annotations]
    base_annotation = annotations[ann_names.index(args.base_annotation)]
    agg_annotations = [annotations[ann_names.index(i)] for i in args.agg_annotations.split(',') if i in ann_names]

    for ann in agg_annotations:
        agged_annotation = spatially_aggregate(ann,[base_annotation],separate=False,summarize=False)
        
        if "/" in ann['properties']['name']:
            ann['properties']['name'] = ann['properties']['name'].replace('/','_')

        export_annotations(
            agged_annotation,
            format='histomics',
            save_path = os.getcwd()+f'/{ann["properties"]["name"]}.json'
        )

        with open(os.getcwd()+f'/{ann["properties"]["name"]}.json','r') as f:
            formatted_anns = json.load(f)
            f.close()

        gc.post(
            f'/annotation/item/{image_id}?token={args.girderToken}',
            data = json.dumps(formatted_anns),
            headers = {
                'X-HTTP-Method':'POST',
                'Content-Type': 'application/json'
            }
        )

    


if __name__=='__main__':
    main(CLIArgumentParser().parse_args())


