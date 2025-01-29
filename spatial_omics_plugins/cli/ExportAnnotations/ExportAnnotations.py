"""Plugin allowing for the exporting of annotations in various formats
"""

import os
import sys
import numpy as np
import json

import girder_client
from ctk_cli import CLIArgumentParser

from fusion_tools.utils.shapes import export_annotations
from fusion_tools.handler.dsa_handler import DSAHandler

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
        girderApiUrl= args.girderApiUrl
    )

    annotations = dsa_handler.get_annotations(
        item = image_id
    )

    export_annotations(
        annotations,
        format = args.output_format,
        save_path = os.getcwd()+f'/{image_name}.{args.output_format}'
    )

    gc.uploadFileToItem(
        itemId = image_id,
        filepath = os.getcwd()+f'/{image_name}.{args.output_format}'
    )


if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

