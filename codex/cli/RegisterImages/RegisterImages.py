"""Plugin for registering two images
"""
import os
import sys
import large_image.constants
import numpy as np
import json
import yaml
import math
import copy


import girder_client
from ctk_cli import CLIArgumentParser

import large_image
import large_image_converter
import pystackreg
import numpy as np
import rasterio.features
import shapely

import tempfile
import logging

import skimage.filters
import skimage.morphology
import skimage.transform
from histomicstk.preprocessing.color_deconvolution import color_deconvolution
from histomicstk.preprocessing.color_deconvolution.stain_color_map import stain_color_map

def annotation_to_shapely(annot, reduce=1):
    return shapely.polygons([
        shapely.linearrings([[p[0] / reduce, p[1] / reduce]
                             for p in element['points']]
        for element in annot['annotation']['elements']
        if element['type']=='polyline' and element['closed'])
    ])

def segment_nuclei_base(image,is_histology):
    
    if is_histology:
        stains = ['hematoxylin','eosin','null']
        stain_matrix = np.array([stain_color_map[stain] for stain in stains]).T

        grayscale_image = color_deconvolution(image, stain_matrix).Stains[:,:,0]
        grayscale_image = np.squeeze(255-grayscale_image)
    
    else:

        if len(image.shape)==3:
            if image.shape[-1]>1:
                grayscale_image = np.mean(image,axis=-1)
            else:
                grayscale_image = np.squeeze(grayscale_image)
        else:
            grayscale_image = image

    thresh_image = (grayscale_image>skimage.filters.threshold_otsu(grayscale_image))
    small_object_removed = skimage.morphology.remove_small_objects(thresh_image,5)
    opened_image = skimage.morphology.binary_opening(small_object_removed,skimage.morphology.disk(3))

    final_mask = 255*skimage.morphology.remove_small_objects(opened_image,5)

    return final_mask

def extract_alignment_base(gc,tile_source, sizeX, sizeY, region, frame, reduce, rescale_val, annotation_args, is_histology):
    """This function extracts the nuclei segmentation regions from an image that will be used in alignment. 
    """

    regionparams = {
        'format': large_image.constants.TILE_FORMAT_NUMPY
    }

    if frame:
        # "Frame" refers to either channels, times, or depth
        # as indicated by IndexC, IndexT, and IndexZ
        regionparams['frame'] = int(frame)

    if region:
        regionparams = regionparams | {'top': region[0], 'left': region[1], 'bottom': region[2], 'right': region[3]}
        sizeX = int(region[3]-region[1])
        sizeY = int(region[2]-region[0])


    # annotation_args will contain information about either segmentation method or annotationID to use
    # if a segmentation method is provided, extract the image region and perform segmentation, otherwise
    # get the annotation from annotationId and create a binary mask
    if 'annotationId' in annotation_args:
        # Pulling that annotationId from a specified region (or just in general)
        if not region:
            nuclei_annotations = gc.get(f'annotation/{annotation_args["annotationId"]}')
        else:
            nuclei_annotations = gc.get(
                f'annotation/{annotation_args["annotationId"]}',
                parameters = {
                    'top': region[0],
                    'left': region[1],
                    'bottom': region[2],
                    'right': region[3]
                }
            )
            sizeX = int(region[3] - region[1])
            sizeY = int(region[2] - region[0])

        shapely_annotations = annotation_to_shapely(nuclei_annotations, reduce / rescale_val)
        
        nuclei_image = (rasterio.features.rasterize(shapely_annotations, out_shape=(sizeY,sizeX))>0).astype('bool')
    else:
        # Segmentation parameters for nuclei segmentation 
        regionparams['output'] = {
            'maxWidth': int(sizeX // reduce*rescale_val),
            'maxHeight': int(sizeY // reduce * rescale_val)
        }

        image = tile_source.getRegion(**regionparams)[0]

        nuclei_image = segment_nuclei_base(image, is_histology)


    return nuclei_image

def register_images(img1,img2,sizeX,sizeY,transform_type,rescale, reduce):
    img1 = np.pad(img1, ((0, sizeY - img1.shape[0]),
                        (0, sizeX - img1.shape[1])), mode='constant')
    sr = pystackreg.StackReg(getattr(
        pystackreg.StackReg, transform_type, pystackreg.StackReg.AFFINE))

    # check four cardinal rotations
    print('Check rotations')
    best = None
    rotMatrix = {
        0: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        1: [[0, -1, img2.shape[1]], [1, 0, 0], [0, 0, 1]],
        2: [[-1, 0, img2.shape[1]], [0, -1, img2.shape[0]], [0, 0, 1]],
        3: [[0, 1, 0], [-1, 0, img2.shape[0]], [0, 0, 1]],
    }
    for rot in range(4):
        rimg = np.rot90(img2, rot)
        rimg = np.pad(rimg, ((0, sizeY - rimg.shape[0]),
                                (0, sizeX - rimg.shape[1])), mode='constant')
        timg = sr.register_transform(img1, rimg)
        print(f'Rotation matrix {rot * 90} deg')
        print(sr.get_matrix())
        if best is None or np.sum(np.abs(timg - img1)) < best[1]:
            best = rot, np.sum(np.abs(timg - img1)), timg, sr.get_matrix().copy()
            print(f'Current best rotation is {rot * 90} deg, score {best[1]}')
    if best[0]:
        print(f'Use rotation of {best[0] * 90} degrees ccw')

    print('Register plain')
    rot = best[0]
    timg = best[2]
    print('Direct result')
    full = best[3]
    print(full)
    print('With rotation')
    full = np.dot(rotMatrix[best[0]], full)
    print(full)
    if rescale != 1:
        full = np.dot([[1. / rescale, 0, 0], [0, 1. / rescale, 0], [0, 0, 1]], full)
    full[0][2] *= reduce
    full[1][2] *= reduce
    print('Full result')
    print(full)
    print('Inverse')
    inv = np.linalg.inv(full)
    print(inv)
    print('Transforming image')

    return full, inv

def transform_images(ts1, ts2, matrix, outmergepath=None, rscale=1):
    
    if hasattr(matrix, 'tolist'):
        matrix = matrix.tolist()
    sx, sy = ts1.sizeX, ts1.sizeY
    for cx, cy in [(0, 0), (ts2.sizeX, 0), (0, ts2.sizeY), (ts2.sizeX, ts2.sizeY)]:
        txy = np.array(matrix).dot([[cx], [cy], [1]])
        sx = int(math.ceil(max(sx, txy[0])))
        sy = int(math.ceil(max(sy, txy[1])))
    trans2 = {
        'name': f'Transform of {os.path.basename(ts2.largeImagePath)}',
        'width': sx,
        'height': sy,
        'backgroundColor': [0, 0, 0],
        'scale': {},
        'dtype': str(ts2.dtype) if np.dtype(ts2.dtype).kind != 'f' else 'uint16',
        'sources': [{
            'path': ts2.largeImagePath,
            'position': {
                's11': matrix[0][0],
                's12': matrix[0][1],
                'x': matrix[0][2],
                's21': matrix[1][0],
                's22': matrix[1][1],
                'y': matrix[1][2],
            },
        }],
    }
    if ts2.frames > 1:
        trans2['singleBand'] = True
    for k in {'mm_x', 'mm_y', 'magnfication'}:
        val = ts2.metadata.get(k) or ts1.metadata.get(k) or 0
        if val:
            trans2['scale'][k] = val / rscale
    print('---')
    print(yaml.dump(trans2, sort_keys=False))

    combo = copy.deepcopy(trans2)
    combo['name'] = (
        f'Transform of {os.path.basename(ts2.largeImagePath)} with '
        f'{os.path.basename(ts1.largeImagePath)}')
    if ts1.frames == 1 and ts2.frames == 1:
        combo['sources'].append({
            'path': ts1.largeImagePath,
            'z': 1,
        })
    elif ts1.frames == 1 and ts2.frames > 1:
        for band in range(1 if ts1.bandCount < 3 else 3):
            combo['sources'].append({
                'path': ts1.largeImagePath,
                'channel': ['red', 'green', 'blue'][band] if ts1.bandCount >= 3 else 'gray',
                'z': 0,
                'c': band + (
                    ts2.frames
                    if ts2.metadata.get('IndexRange', {}).get('IndexC', 0) <= 1 else
                    ts2.metadata['IndexRange']['IndexC']),
                'style': {'dtype': 'uint8', 'bands': [{'band': band + 1, 'palette': 'white'}]},
            })
    else:
        combo['singleBand'] = True
        if ts2.frames == 1:
            src = combo['sources'].pop()
            for band in range(1 if ts2.bandCount < 3 else 3):
                src = src.copy()
                src.update({
                    'channel': ['red', 'green', 'blue'][band] if ts2.bandCount >= 3 else 'gray',
                    'z': 0,
                    'c': band,
                    'style': {'dtype': 'uint8', 'bands': [{'band': band + 1, 'palette': 'white'}]},
                })
                combo['sources'].append(src)
            if 'channels' in ts1.metadata:
                combo['channels'] = [
                    s['channel'] for s in combo['sources']] + ts1.metadata['channels']
        combo['sources'].append({
            'path': ts1.largeImagePath,
            'z': 0,
            'c': (len(combo['sources'])
                  if ts2.metadata.get('IndexRange', {}).get('IndexC', 0) <= 1 else
                  ts2.metadata['IndexRange']['IndexC']),
        })
    print('---')
    print(yaml.dump(combo, sort_keys=False))
    print('\n---')
    with tempfile.TemporaryDirectory() as tmpdir:
        sys.stdout.flush()
        logger = logging.getLogger('large_image')
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        logger.addHandler(handler)
        logger = logging.getLogger('large-image-converter')
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        logger.addHandler(handler)

        if outmergepath:
            output_image(tmpdir, 'outmergetransform.yaml', combo, outmergepath)

def output_image(tmpdir, name, multi, outpath):
    if outpath.rsplit('.')[-1].lower() in {'yaml', 'yml'}:
        multi = copy.deepcopy(multi)
        for source in multi['sources']:
            source['path'] = os.path.join('.', source['path'].rsplit(os.path.sep, 1)[-1])
        open(outpath, 'w').write(yaml.dump(multi))
    else:
        yamlpath = os.path.join(tmpdir, name)
        open(yamlpath, 'w').write(yaml.dump(multi))
        large_image_converter.convert(yamlpath, outpath)

def merge_dict(a:dict, b:dict, path = []):

    # https://stackoverflow.com/questions/7204805/deep-merge-dictionaries-of-dictionaries-in-python
    for key in b:
        if key in a:
            if isinstance(a[key],dict) and isinstance(b[key],dict):
                merge_dict(a[key],b[key],path+[str(key)])
            elif a[key] != b[key]:
                raise Exception(f'Conflict at {".".join(path+[str(key)])}')
        else:
            a[key] = b[key]
    return a

def add_alignment_metadata(gc, item_info, add_meta):

    item_meta = merge_dict(item_info['meta'],add_meta)

    gc.put(
        f'/item/{item_info["_id"]}/metadata',
        parameters = {
            'metadata': item_meta
        }
    )



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

    # Grab both images and move them to the current working directory
    base_file_info = gc.get(f'/file/{args.base_file}')
    trans_file_info = gc.get(f'/file/{args.trans_file}')

    base_file_name = base_file_info['name']
    trans_file_name = trans_file_info['name']

    gc.downloadFile(
        fileId = args.base_file,
        path = os.getcwd()+f'/{base_file_name}'
    )
    gc.downloadFile(
        fileId = args.trans_image,
        path = os.getcwd()+f'/{trans_file_name}'
    )

    base_image = large_image.open(os.getcwd()+f'/{base_file_name}')
    trans_image = large_image.open(os.getcwd()+f'/{trans_file_name}')

    # Reading in downsampled base_image and trans_image based on specifications
    #   # if specific frame in either is meant to be used (multi-channel)
    #   # if specific z-index in either is meant to be used (multi-level)
    #   # if specific region is specified for either

    # Original rescale value
    rescale_val = 1.0
    # Whether or not the microns per-pixel (MPP) value is present in both images (corrects for varying scales)
    mpp_present = base_image.metadata['mm_x'] and base_image.metadata['mm_x'] and trans_image.metadata['mm_x'] and trans_image['mm_y']
    if args.trans_image_scale and mpp_present:
        # Calculate initial rescaling based on differences in MPP
        rescale_val = math.sqrt(
            trans_image.metadata['mm_x'] * trans_image.metadata['mm_y'] /
            (base_image.metadata['mm_x'] * base_image.metadata['mm_y'])
        )
    maxRes = int(max(base_image.sizeX, base_image.sizeY, trans_image.sizeX * rescale_val, trans_image.sizeY * rescale_val))
    reduce = 1
    while math.ceil(maxRes / reduce) >= args.downsample_dim:
        reduce *= 2
    
    sizeX = int(math.ceil(max(base_image.sizeX, trans_image.sizeX * rescale_val, trans_image.sizeY * rescale_val) / reduce))
    sizeY = int(math.ceil(max(base_image.sizeY, trans_image.sizeX * rescale_val, trans_image.sizeY * rescale_val) / reduce))

    # Read in cell annotations if they are provided

    # If not provided, segment cells in downsampled images
    #   # specify image type (histology, IF (multiframe))
    #   # segmentation method based on type or specification
    base_image_alignment_image = extract_alignment_base(
        gc = gc,
        tile_source = base_image,
        sizeX = sizeX,
        sizeY = sizeY,
        region = args.region1 if args.region1 else None,
        frame = args.frame_index_1 if args.use_frame_index_1 else None,
        reduce = reduce,
        rescale_val = rescale_val,
        annotation_args = {
            'annotationId': args.annotation_id_1
        } if not args.annotation_id_1 is None else {}
    )

    trans_image_alignment_image = extract_alignment_base(
        gc = gc,
        tile_source = trans_image,
        sizeX = sizeX,
        sizeY = sizeY,
        region = args.region2 if args.region2 else None,
        frame = args.frame_index_1 if args.use_frame_index_2 else None,
        reduce = reduce,
        rescale_val = rescale_val,
        annotation_args = {
            'annotationId': args.annotation_id_2
        } if not args.annotation_id_2 is None else {}
    )

    full_matrix, inverse_matrix = register_images(
        base_image_alignment_image,
        trans_image_alignment_image,
        sizeX,
        sizeY,
        args.transform
    )

    # posting the matrices to item metadata
    base_meta = {
        'Alignment': {
            'inverse': inverse_matrix.astype(float).tolist(),
            'full': full_matrix.astype(float).tolist(),
            '_id': trans_file_info['itemId'],
            'name': trans_file_info['name']
        }
    }

    base_item_info = gc.get(f'/item/{base_file_info["itemId"]}')
    add_alignment_metadata(gc,base_item_info,base_meta)

    # Creating the merged image if specified
    if args.merge_images:
        transform_images(
            base_image,
            trans_image,
            inverse_matrix,
            outmergepath=args.output_merged_image,
            rscale = rescale_val
        )

if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

