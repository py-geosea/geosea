# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
#  Purpose: Convenience imports for geosea
#   Author: Florian Petersen
#           Katrin Hannemann
#
# GEOMAR Helmholtz Centre for Ocean Research Kiel, Germany
#
# Version: 1.21                 August 2020
#
# -----------------------------------------------------------------------------
"""
GeoSEA: A Python Toolbox for seafloor geodesy
======================================================

GeoSEA is an open-source project to provide a pursuning Python tool for
seafloor gedetic data processing.
"""

from __future__ import absolute_import
import warnings

from .read import *
from .read_id import *
from .read_data import *
from .read_bsl import *
from .read_meta import *
from .read_tides import *
from .read_airpressure import *

from .proc_bsl import *

from .extract_df import *
from .search_df import *
from .create_df import *

from .range_sv import *
from .sw import *
from .replace import *


from .vert_bsl import *
from .hori_bsl import *

from .change2dateindex import *
from .change_dtype import *

from .compare_df import *

from .calc import *

global GMT_DATEFORMAT # Output date format
global IN_DATEFORMAT # Input date format
global PROJECTS # GeoSEA projects

GMT_DATEFORMAT = '%Y-%m-%dT%H:%M'
IN_DATEFORMAT = '%Y/%m/%d %H:%M:%S'

# MAR = MARSITE
# CHI = GeoSEA
# ETN = MARGOMET
PROJECTS = {'MAR' : '2014-11-16 00:00:00', 'CHI' : '2015-12-14 00:00:00', 'ETN' : '2016-04-15 00:00:00'}

