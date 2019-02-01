import os
import sys
import ssim
import shutil

def findsplat(workdir):
    """Starting at workdir, search each successive parent dir until we find
       the splat executable"""
    (workdir, tail) = os.path.split(workdir)
    while tail:
        (workdir, tail) = os.path.split(workdir)
        path = shutil.which("splat", path=workdir)
        if path:
            return path

def pyssim(base_image, comparison_image):
    """Starting at workdir, search each successive parent dir until we find
       the splat executable"""
    gaussian_kernel_sigma = 1.5
    gaussian_kernel_width = 11
    gaussian_kernel_1d = ssim.get_gaussian_kernel(
            gaussian_kernel_width, gaussian_kernel_sigma)

    base_ss = ssim.SSIM(base_image, gaussian_kernel_1d)
    ssim_value = base_ss.ssim_value(comparison_image)
    #sys.stdout.write('%.7g' % ssim_value);
    return ssim_value

def striplinefromfile(filename, text):
    f = open(filename,"r+",encoding = "ISO-8859-1")
    data = f.readlines()
    f.seek(0)
    for i in data:
        if not text in i:
            f.write(i)
    f.truncate()
    f.close()
