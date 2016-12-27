import os.path as osp
import SimpleITK as sitk

from .dicom import *

def image_read(filename, **args):
    if osp.isdir(filename):
        image = load_dicom(filename)
    else:
        image = sitk.ReadImage(filename)
    return image

def image_write(image, filename, **args):
    use_compression = True
    sitk.WriteImage(image, filename, use_compression)
    return
