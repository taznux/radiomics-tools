import os.path as osp

pkg_dir = osp.abspath(osp.dirname(__file__))

from . import io
from . import util
from .io import metadata as metadata
