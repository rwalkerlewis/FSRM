# Velocity Models Directory
#
# 1-D layered velocity models for seismic wave propagation in FSRM.
#
# Models are defined as Python classes in scripts/fsrm/velocity_models.py.
# This directory holds supplementary data files (e.g. depth profiles from
# published studies, TauP .npz format, or raw digitised curves) that can
# be loaded by the Python or C++ code.
#
# Available models:
#   punggye_ri  — Korean Peninsula (DPRK test site, Moho 33 km)
#   lop_nor     — Kuruktag / Tarim Basin (China, Moho 48 km)
#   nts         — Nevada Test Site (USA, Moho 30 km)
#   generic     — Generic continental granite (Moho 35 km)
#
# To add a new model, subclass VelocityModel1D in velocity_models.py and
# register it in the _MODELS dict.
