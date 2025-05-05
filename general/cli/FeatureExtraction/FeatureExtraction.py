"""Plugin for programmatic feature extraction
"""
import os
import sys
import numpy as np
import json
import pandas as pd
import requests
import time

import girder_client
from ctk_cli import CLIArgumentParser

from fusion_tools.feature_extraction import ParallelFeatureExtractor

import geopandas as gpd
from shapely.geometry import shape
from fusion_tools.handler.dsa_handler import DSAHandler
from fusion_tools.utils.shapes import histomics_to_geojson

from skimage.segmentation import watershed
from skimage import exposure
from skimage.color import rgb2hsv
import scipy.ndimage as ndi

from skimage.morphology import remove_small_holes, remove_small_objects
from skimage.feature import peak_local_max


def stain_mask(image,mask, seg_params=None):

    if seg_params is None:
        seg_params = [
            {
                'name': 'Nuclei',
                'threshold': 150,
                'min_size': 40
            },
            {
                'name': 'Eosinophilic',
                'threshold': 30,
                'min_size': 20
            },
            {
                'name': 'Luminal Space',
                'threshold': 0,
                'min_size': 0
            }
        ]

    image_shape = np.shape(image)

    sub_comp_image = np.zeros((image_shape[0],image_shape[1],3))
    remainder_mask = np.ones((image_shape[0],image_shape[1]))

    hsv_image = np.uint8(255*rgb2hsv(image))
    hsv_image = hsv_image[:,:,1]

    for idx,param in enumerate(seg_params):

        # Check for if the current sub-compartment is nuclei
        if param['name'].lower()=='nuclei':
            # Using the inverse of the value channel for nuclei
            h_image = 255-np.uint8(255*rgb2hsv(image)[:,:,2])
            h_image = np.uint8(255*exposure.equalize_hist(h_image))

            remaining_pixels = np.multiply(h_image,remainder_mask)
            masked_remaining_pixels = np.multiply(remaining_pixels,mask)
            #masked_remaining_pixels = remaining_pixels

            # Applying manual threshold
            masked_remaining_pixels[masked_remaining_pixels<=param['threshold']] = 0
            masked_remaining_pixels[masked_remaining_pixels>0] = 1

            # Area threshold for holes is controllable for this
            sub_mask = remove_small_holes(masked_remaining_pixels>0,area_threshold=10)
            sub_mask = sub_mask>0
            # Watershed implementation from: https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_watershed.html
            distance = ndi.distance_transform_edt(sub_mask)
            labeled_mask, _ = ndi.label(sub_mask)
            coords = peak_local_max(distance,footprint=np.ones((3,3)),labels = labeled_mask)
            watershed_mask = np.zeros(distance.shape,dtype=bool)
            watershed_mask[tuple(coords.T)] = True
            markers, _ = ndi.label(watershed_mask)
            sub_mask = watershed(-distance,markers,mask=sub_mask)
            sub_mask = sub_mask>0

            # Filtering out small objects again
            sub_mask = remove_small_objects(sub_mask,param['min_size'])

        else:

            remaining_pixels = np.multiply(hsv_image,remainder_mask)
            masked_remaining_pixels = np.multiply(remaining_pixels,mask)
            #masked_remaining_pixels = remaining_pixels

            # Applying manual threshold
            masked_remaining_pixels[masked_remaining_pixels<=param['threshold']] = 0
            masked_remaining_pixels[masked_remaining_pixels>0] = 1

            # Filtering by minimum size
            small_object_filtered = (1/255)*np.uint8(remove_small_objects(masked_remaining_pixels>0,param['min_size']))

            sub_mask = small_object_filtered

        sub_comp_image[sub_mask>0,idx] = 1
        remainder_mask -= sub_mask>0

    # Assigning remaining pixels within the boundary mask to the last sub-compartment
    remaining_pixels = np.multiply(mask,remainder_mask)
    #remaining_pixels = remainder_mask
    sub_comp_image[remaining_pixels>0,idx] = 1

    final_mask = np.zeros_like(remainder_mask)
    final_mask += sub_comp_image[:,:,0]
    final_mask += 2*sub_comp_image[:,:,1]
    final_mask += 3*sub_comp_image[:,:,2]

    return final_mask



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

    dsa_handler = DSAHandler(
        girderApiUrl=args.girderApiUrl
    )

    tile_server = dsa_handler.get_tile_server(
        item = image_id
    )

    annotations = []
    for attempt in range(10):
        print(f'On attempt: {attempt}/10')
        try:
            annotations = dsa_handler.get_annotations(
                item = image_id,
                user_token = args.girderToken
            )
            break
        except Exception as e:
            print(e)
            print('-----------------')
            print('Retrying')
            print('-----------------')
            # This will either be a requests.exceptions.ChunkedEncodingError or urllib3.exceptions.ProtocolError or both 
            # but not sure if this will also trigger a girder_client.HttpError
            time.sleep(1)

    # Backup method for getting annotations
    if len(annotations)==0:
        print('Loading annotations another way')
        histomics_annotations = dsa_handler.gc.get(f'annotation/item/{image_id}')
        annotations = histomics_to_geojson(histomics_annotations)


    if not args.extract_sub_compartments:

        feature_extractor = ParallelFeatureExtractor(
            image_source = tile_server,
            feature_list = [
                'distance_transform',
                'morphology',
                'color',
                'texture'
            ],
            preprocess = None,
            n_jobs = 4,
            verbose = True
        )
    else:

        #TODO: Have to define these parameters somewhere
        mask_names = ["Nuclei","Eosinophilic","Luminal Space"]
        channel_names = ["Red","Green","Blue"]

        sub_comp_params = [
            {
                'name': 'Nuclei',
                'threshold': args.hematoxylin_threshold,
                'min_size': args.hematoxylin_min_size
            },
            {
                'name': 'Eosinophilic',
                'threshold': args.eosinophilic_threshold,
                'min_size': args.eosinophilic_min_size
            },
            {
                'name': 'Luminal Space',
                'threshold': 0,
                'min_size': 0
            }
        ]


        feature_extractor = ParallelFeatureExtractor(
            image_source = tile_server,
            feature_list = [
                'distance_transform',
                'morphology',
                'color',
                'texture'
            ],
            preprocess = None,
            sub_mask = lambda image,mask: stain_mask(image,mask,seg_params = sub_comp_params),
            mask_names = mask_names,
            channel_names = channel_names,
            n_jobs = 4,
            verbose = True
        )

    # Add these dataframes as files? Or store in the "user" field?
    histomics_anns = []
    for ann_idx,ann in enumerate(annotations):
        print(f'On annotation: {ann_idx+1}/{len(annotations)}')
        feature_df = feature_extractor.start(ann['features'])
        feature_df = pd.json_normalize(feature_df.to_dict('records')).dropna(how='all').fillna(0)
        histomics_ann = {
            'annotation': {
                'name': ann['properties']['name'],
                '_id': ann['properties']['_id'],
                'elements': []
            }
        }

        # Feature_df will have bounding box columns (min_x, min_y, max_x, max_y) which can be used to align with structures
        if args.save_to_files:
            save_path = os.getcwd()+f'/{ann["properties"]["name"].replace("/","_").replace(".","_")}.csv'
            feature_df.to_csv(save_path)
            gc.uploadFileToItem(image_id, save_path, reference=None, mimeType=None, filename=None, progressCallback=None)
        
        if args.save_to_elements:

            if not feature_df.empty:
                # They should just match by index but this will double check to ensure that bounding boxes align between Features (in the geojson) and feature rows in the dataframe
                bbox_coords = feature_df[['bbox.min_x','bbox.min_y','bbox.max_x','bbox.max_y']].values.tolist()

                for f in ann['features']:
                    ann_bbox = list(shape(f['geometry']).bounds)
                    feature_df_idx = bbox_coords.index(ann_bbox)
                    f['properties'] = f['properties'] | feature_df.iloc[feature_df_idx,:].to_dict()
                    histomics_el = {
                        'type': 'polyline',
                        'points': [i+[0] if len(i)==2 else i for i in f['geometry']['coordinates'][0]],
                        'user': f['properties']
                    }

                    histomics_ann['annotation']['elements'].append(histomics_el)

        histomics_anns.append(histomics_ann)
            
    if args.save_to_elements:
        for h in histomics_anns:
            put_dict = {
                "name": h['annotation']['name'],
                "elements": h['annotation']['elements']
            }
            gc.put(
                f'/annotation/{h["annotation"]["_id"]}?token={args.girderToken}',
                data = json.dumps(put_dict)
            )

    # Putting job parameters to item metadata:
    job_submitter = gc.get('/user/me')
    job_meta = {
        'input_image': args.input_image,
        'extract_sub_compartments': args.extract_sub_compartments,
        'save_to_elements': args.save_to_elements,
        'save_to_files': args.save_to_files,
        'user': job_submitter['login']
    }

    if args.extract_sub_compartments:
        job_meta = job_meta | {
            'eosinophilic_min_size': args.eosinophilic_min_size,
            'eosinophilic_threshold': args.eosinophilic_threshold,
            'hematoxylin_min_size': args.hematoxylin_min_size,
            'hematoxylin_threshold': args.hematoxylin_threshold
        }

    gc.put(f'/item/{image_id}/metadata',parameters={'metadata': json.dumps(job_meta)})

def test():
    
    test_args = {
        'input_image': '6717e743f433060d2884838c',
        'girderApiUrl': 'http://ec2-3-230-122-132.compute-1.amazonaws.com:8080/api/v1/',
        'girderToken': 'MG87Dcak28A218XA4HxBLOd1KRaF3s8Vq6YHf79sbQqnaL0hBwPP3aYRsTObQ6tZ',
        'eosinophilic_min_size': 20,
        'eosinophilic_threshold': 30,
        'extract_sub_compartments': True,
        'hematoxylin_min_size': 40,
        'hematoxylin_threshold': 150,
        'save_to_elements': True,
        'save_to_files': True
    }
    class ArgsObj:
        def __init__(self, args_dict):
            for a in args_dict:
                setattr(self,a,args_dict[a])

    main(ArgsObj(test_args))



if __name__=='__main__':
    main(CLIArgumentParser().parse_args())
    #test()
