import SimpleITK as sitk

def image_read(filename, **args):
    image = sitk.ReadImage(filename)
    return image

def image_write(image, filename, **args):
    use_compression = True
    sitk.WriteImage(image, filename, use_compression)
    return
