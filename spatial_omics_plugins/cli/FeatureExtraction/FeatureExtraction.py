"""Plugin for programmatic feature extraction
"""
import os
import sys
import numpy as np
import json

import girder_client
from ctk_cli import CLIArgumentParser

from fusion_tools.feature_extraction import (
    ParallelFeatureExtractor, distance_transform_features, color_features,
    morphological_features, texture_features
)

from fusion_tools.tileserver import DSATileServer
from fusion_tools.handler import DSAHandler

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
    
    image_id = gc.get(f'/file/{args.input_image}')['itemId']

    tile_server = DSATileServer(
        api_url = args.girderApiUrl,
        item_id = image_id
    )

    annotations = DSAHandler(
        girderApiUrl=args.girderApiUrl
    ).get_annotations(
        item = image_id
    )

    feature_extractor = ParallelFeatureExtractor(
        image_source = tile_server,
        feature_list = [
            lambda image, mask, coords: distance_transform_features(image, mask, coords),
            lambda image, mask, coords: color_features(image, mask, coords),
            lambda image, mask, coords: texture_features(image, mask, coords),
            lambda image, mask, coords: morphological_features(image, mask, coords)
        ],
        preprocess = None,
        n_jobs = 4,
        verbose = True
    )

    for ann in annotations:
        feature_df = feature_extractor.start(ann['features'])



if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

