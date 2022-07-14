# -*- coding: utf-8 -*-
# Copyright 2018-2022 the orix developers
#
# This file is part of orix.
#
# orix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# orix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with orix.  If not, see <http://www.gnu.org/licenses/>.

"""Reader of a crystal map from an .ang text file in formats produced by
EDAX TSL, NanoMegas ASTAR Index, or EMsoft's EMdpmerge program.
"""

import re
import warnings

from diffpy.structure import Lattice, Structure
import numpy as np

from orix.crystal_map import CrystalMap, PhaseList, create_coordinate_arrays
from orix.quaternion.rotation import Rotation
from orix.quaternion.symmetry import point_group_aliases
from orix import __version__


# MTEX has this format sorted out, check out their readers when fixing
# issues and adapting to other versions of this file format in the future:
# https://github.com/mtex-toolbox/mtex/blob/develop/interfaces/loadEBSD_ang.m
# https://github.com/mtex-toolbox/mtex/blob/develop/interfaces/loadEBSD_ACOM.m

# Plugin description
format_name = "ang"
file_extensions = ["ang"]
writes = True
writes_this = CrystalMap


def file_reader(filename):
    """Return a :class:`~orix.crystal_map.crystal_map.CrystalMap` object
    from a file in EDAX TLS's .ang format. The map in the input is
    assumed to be 2D.

    Many vendors produce an .ang file. Supported vendors are
    * EDAX TSL

    * NanoMegas ASTAR Index

    * EMsoft (from program `EMdpmerge`)

    * orix

    All points satisfying the following criteria are classified as not
    indexed:

    * EDAX TSL: confidence index == -1

    Parameters
    ----------
    filename : str
        Path and file name.

    Returns
    -------
    CrystalMap
    """
    # Get file header
    with open(filename) as f:
        header = _get_header(f)

    # Get phase names and crystal symmetries from header (potentially empty)
    phase_ids, phase_names, symmetries, lattice_constants = _get_phases_from_header(
        header
    )
    structures = []
    for name, abcABG in zip(phase_names, lattice_constants):
        structures.append(Structure(title=name, lattice=Lattice(*abcABG)))

    # Read all file data
    file_data = np.loadtxt(filename)

    # Get vendor and column names
    n_rows, n_cols = file_data.shape
    vendor, column_names = _get_vendor_columns(header, n_cols)

    # Data needed to create a CrystalMap object
    data_dict = {
        "euler1": None,
        "euler2": None,
        "euler3": None,
        "x": None,
        "y": None,
        "phase_id": None,
        "prop": {},
    }
    for column, name in enumerate(column_names):
        if name in data_dict.keys():
            data_dict[name] = file_data[:, column]
        else:
            data_dict["prop"][name] = file_data[:, column]

    # Add phase list to dictionary
    data_dict["phase_list"] = PhaseList(
        names=phase_names,
        point_groups=symmetries,
        structures=structures,
        ids=phase_ids,
    )

    # Set which data points are not indexed
    if vendor in ["orix", "tsl"]:
        data_dict["phase_id"][np.where(data_dict["prop"]["ci"] == -1)] = -1
    # TODO: Add not-indexed convention for INDEX ASTAR

    # Set scan unit
    if vendor in ["tsl", "emsoft"]:
        scan_unit = "um"
    else:  # NanoMegas
        scan_unit = "nm"
    data_dict["scan_unit"] = scan_unit

    # Create rotations
    data_dict["rotations"] = Rotation.from_euler(
        np.column_stack(
            (data_dict.pop("euler1"), data_dict.pop("euler2"), data_dict.pop("euler3"))
        ),
    )

    return CrystalMap(**data_dict)


def _get_header(file):
    """Return the first lines starting with '#' in an .ang file.

    Parameters
    ----------
    file : _io.TextIO
        File object.

    Returns
    -------
    header : list
        List with header lines as individual elements.
    """
    header = []
    line = file.readline()
    while line.startswith("#"):
        header.append(line.rstrip())
        line = file.readline()
    return header


def _get_vendor_columns(header, n_cols_file):
    """Return the .ang file column names and vendor, determined from the
    header.

    Parameters
    ----------
    header : list
        List with header lines as individual elements.
    n_cols_file : int
        Number of file columns.

    Returns
    -------
    vendor : str
        Determined vendor ("tsl", "astar", or "emsoft").
    column_names : list of str
        List of column names.
    """
    # Assume EDAX TSL by default
    vendor = "tsl"

    # Determine vendor by searching for the vendor footprint in the header
    vendor_footprint = {
        "emsoft": "EMsoft",
        "astar": "ACOM",
        "orix": "Column names: phi1, Phi, phi2",
    }
    footprint_line = None
    for name, footprint in vendor_footprint.items():
        for line in header:
            if footprint in line:
                vendor = name
                footprint_line = line
                break

    # Vendor column names
    column_names = {
        "unknown": [
            "euler1",
            "euler2",
            "euler3",
            "x",
            "y",
            "unknown1",
            "unknown2",
            "phase_id",
        ],
        "tsl": [
            "euler1",
            "euler2",
            "euler3",
            "x",
            "y",
            "iq",  # Image quality from Hough transform
            "ci",  # Confidence index
            "phase_id",
            "unknown1",
            "fit",  # Pattern fit
            "unknown2",
            "unknown3",
            "unknown4",
            "unknown5",
        ],
        "emsoft": [
            "euler1",
            "euler2",
            "euler3",
            "x",
            "y",
            "iq",  # Image quality from Krieger Lassen's method
            "dp",  # Dot product
            "phase_id",
        ],
        "astar": [
            "euler1",
            "euler2",
            "euler3",
            "x",
            "y",
            "ind",  # Correlation index
            "rel",  # Reliability
            "phase_id",
            "relx100",  # Reliability x 100
        ],
        "orix": [
            "euler1",
            "euler2",
            "euler3",
            "x",
            "y",
            "iq",
            "ci",
            "phase_id",
            "detector_signal",
            "fit",
        ],
    }

    n_cols_expected = len(column_names[vendor])
    if vendor == "orix" and "Column names" in footprint_line:
        # Append names of extra properties found, if any, in the orix
        # .ang file header
        n_cols = len(column_names[vendor])
        extra_props = footprint_line.split(":")[1].split(",")[n_cols:]
        column_names[vendor] += [i.lstrip(" ").replace(" ", "_") for i in extra_props]
    elif n_cols_file != n_cols_expected:
        warnings.warn(
            f"Number of columns, {n_cols_file}, in the file is not equal to "
            f"the expected number of columns, {n_cols_expected}, for the \n"
            f"assumed vendor '{vendor}'. Will therefore assume the following "
            "columns: euler1, euler2, euler3, x, y, unknown1, unknown2, "
            "phase_id, unknown3, unknown4, etc."
        )
        vendor = "unknown"
        n_cols_unknown = len(column_names["unknown"])
        if n_cols_file > n_cols_unknown:
            # Add potential extra columns to properties
            for i in range(n_cols_file - n_cols_unknown):
                column_names["unknown"].append("unknown" + str(i + 3))

    return vendor, column_names[vendor]


def _get_phases_from_header(header):
    """Return phase names and symmetries detected in an .ang file
    header.

    Parameters
    ----------
    header : list
        List with header lines as individual elements.

    Returns
    -------
    ids : list of int
        Phase IDs.
    phase_names : list of str
        List of names of detected phases.
    phase_point_groups : list of str
        List of point groups of detected phase.
    lattice_constants : list of list of floats
        List of list of lattice parameters of detected phases.

    Notes
    -----
    Regular expressions are used to collect phase name, formula and
    point group. This function have been tested with files from the
    following vendor's formats: EDAX TSL OIM Data Collection v7, ASTAR
    Index, and EMsoft v4/v5.
    """
    regexps = {
        "id": "# Phase([ \t]+)([0-9 ]+)",
        "name": "# MaterialName([ \t]+)([A-z0-9 ]+)",
        "formula": "# Formula([ \t]+)([A-z0-9 ]+)",
        "point_group": "# Symmetry([ \t]+)([A-z0-9 ]+)",
        "lattice_constants": r"# LatticeConstants([ \t+])(.*)",
    }
    phases = {
        "name": [],
        "formula": [],
        "point_group": [],
        "lattice_constants": [],
        "id": [],
    }
    for line in header:
        for key, exp in regexps.items():
            match = re.search(exp, line)
            if match:
                group = re.split("[ \t]", line.lstrip("# ").rstrip(" "))
                group = list(filter(None, group))
                if key == "name":
                    group = " ".join(group[1:])  # Drop "MaterialName"
                elif key == "lattice_constants":
                    group = [float(i) for i in group[1:]]
                else:
                    group = group[-1]
                phases[key].append(group)

    # Check if formula is empty (sometimes the case for ASTAR Index)
    names = phases["formula"]
    if len(names) == 0 or any([i != "" for i in names]):
        names = phases["name"]

    # Ensure each phase has an ID (hopefully found in the header)
    phase_ids = [int(i) for i in phases["id"]]
    n_phases = len(phases["name"])
    if len(phase_ids) == 0:
        phase_ids += [i for i in range(n_phases)]
    elif n_phases - len(phase_ids) > 0 and len(phase_ids) != 0:
        next_id = max(phase_ids) + 1
        n_left = n_phases - len(phase_ids)
        phase_ids += [i for i in range(next_id, next_id + n_left)]

    # Check header symmetry notation then convert and return corresponding laue groups
    phases["point_group"] = _symmetry_unifier(phases)

    return phase_ids, names, phases["point_group"], phases["lattice_constants"]


def _symmetry_unifier(phases):
    """ Check symmetry notation uniformity and return laue groups
    Parameters
    ----------
    phases : dict
        dict of phase information extracted from .ang header
    Returns
    -------
    phases["point_group"] : list
        laue group of phases
    """
    # Dictionary of accepted International, TSL, Schoenflies, and Geology notation point group names
    symmetry_notations = ["International notation", "TSL", "Schoenflies", "Geology"]
    # laue_group_dict[notation][symmetry] returns laue group for given notations crystal symmetry
    laue_group_dict = {
        "International notation": {'1': '-1', '-1': '-1', '2': '2/m', 'm': '2/m', '2/m': '2/m', '222': 'mmm',
                                   'mm2': 'mmm', 'mmm': 'mmm', '4': '4/m', '-4': '4/m', '4/m': '4/m', '422': '4/mmm',
                                   '4mm': '4/mmm', '42m': '4/mmm', '4/mmm': '4/mmm', '3': '-3', '-3': '-3', '32': '-3m',
                                   '3m': '-3m', '-3m': '-3m', '6': '6/m', '-6': '6/m', '6/m': '6/m', '622': '6/mmm',
                                   '6mm': '6/mmm', '-6m2': '6/mmm', '6/mmm': '6/mmm', '23': 'm-3', 'm-3': 'm-3',
                                   '432': 'm-3m', '-43m': 'm-3m', 'm-3m': 'm-3m'},
        "TSL": {'1': '-1', '20': '2/m', '22': 'mmm', '4': '4/m', '42': '4/mmm', '3': '-3', '32': '-3m', '6': '6/m',
                '62': '6/mmm', '23': 'm-3', '43': 'm-3m'},
        "Schoenflies": {'C1': '-1', 'S2': '-1', 'C2': '2/m', 'C1h': '2/m', 'C2h': '2/m', 'V': 'mmm', 'C2v': 'mmm',
                        'D2h': 'mmm', 'C4': '4/m', 'S4': '4/m', 'C4h': '4/m', 'D4': '4/mmm', 'C4v': '4/mmm',
                        'D2d': '4/mmm', 'D4h': '4/mmm', 'C3': '-3', 'S6': '-3', 'D3': '-3m', 'C3v': '-3m', 'D3d': '-3m',
                        'C6': '6/m', 'C3h': '6/m', 'C6h': '6/m', 'D6': '6/mmm', 'C6v': '6/mmm', 'D3h': '6/mmm',
                        'D6h': '6/mmm', 'T': 'm-3', 'Th': 'm-3', 'O': 'm-3m', 'Td': 'm-3m', 'Oh': 'm-3m'},
        "Geology": {'-1': '-1', '-(22)': '-1', '-2': '2/m', '1': '2/m', '-22': '2/m', '-2-2': 'mmm', '2': 'mmm',
                    '22': 'mmm', '-4': '4/m', '-(42)': '4/m', '-42': '4/m', '-4-2': '4/mmm', '4': '4/mmm',
                    '4-2': '4/mmm', '42': '4/mmm', '-3': '-3', '-(62)': '-3', '-3-2': '-3m', '3': '-3m', '6-2': '-3m',
                    '-6': '6/m', '-32': '6/m', '-62': '6/m', '-6-2': '6/mmm', '6': '6/mmm', '32': '6/mmm',
                    '62': '6/mmm', '-3-3': 'm-3', '4-3': 'm-3', '-4-3': 'm-3m', '-33': 'm-3m', '-43': 'm-3m'},
        "Point Group": {'1': '-1', '2': '-1', '3': '2/m', '4': '2/m', '5': '2/m', '6': 'mmm',
                                   '7': 'mmm', '8': 'mmm', '9': '4/m', '10': '4/m', '11': '4/m', '12': '4/mmm',
                                   '13': '4/mmm', '14': '4/mmm', '15': '4/mmm', '16': '-3', '17': '-3', '18': '-3m',
                                   '19': '-3m', '20': '-3m', '21': '6/m', '22': '6/m', '23': '6/m', '24': '6/mmm',
                                   '25': '6/mmm', '26': '6/mmm', '27': '6/mmm', '28': 'm-3', '29': 'm-3',
                                   '30': 'm-3m', '31': 'm-3m', '32': 'm-3m'},
        "Space Group": {'1': '-1', '2': '-1', '3': '-1', '4': '2/m', '5': '2/m', '6': '2/m', '7': '2/m', '8': '2/m',
                        '9': '2/m', '10': '2/m', '11': '2/m', '12': '2/m', '13': '2/m', '14': '2/m', '15': '2/m',
                        '16': '2/m', '17': 'mmm', '18': 'mmm', '19': 'mmm', '20': 'mmm', '21': 'mmm', '22': 'mmm',
                        '23': 'mmm', '24': 'mmm', '25': 'mmm', '26': 'mmm', '27': 'mmm', '28': 'mmm', '29': 'mmm',
                        '30': 'mmm', '31': 'mmm', '32': 'mmm', '33': 'mmm', '34': 'mmm', '35': 'mmm', '36': 'mmm',
                        '37': 'mmm', '38': 'mmm', '39': 'mmm', '40': 'mmm', '41': 'mmm', '42': 'mmm', '43': 'mmm',
                        '44': 'mmm', '45': 'mmm', '46': 'mmm', '47': 'mmm', '48': 'mmm', '49': 'mmm', '50': 'mmm',
                        '51': 'mmm', '52': 'mmm', '53': 'mmm', '54': 'mmm', '55': 'mmm', '56': 'mmm', '57': 'mmm',
                        '58': 'mmm', '59': 'mmm', '60': 'mmm', '61': 'mmm', '62': 'mmm', '63': 'mmm', '64': 'mmm',
                        '65': 'mmm', '66': 'mmm', '67': 'mmm', '68': 'mmm', '69': 'mmm', '70': 'mmm', '71': 'mmm',
                        '72': 'mmm', '73': 'mmm', '74': 'mmm', '75': 'mmm', '76': '4/m', '77': '4/m', '78': '4/m',
                        '79': '4/m', '80': '4/m', '81': '4/m', '82': '4/m', '83': '4/m', '84': '4/m', '85': '4/m',
                        '86': '4/m', '87': '4/m', '88': '4/m', '89': '4/m', '90': '4/mmm', '91': '4/mmm', '92': '4/mmm',
                        '93': '4/mmm', '94': '4/mmm', '95': '4/mmm', '96': '4/mmm', '97': '4/mmm', '98': '4/mmm',
                        '99': '4/mmm', '100': '4/mmm', '101': '4/mmm', '102': '4/mmm', '103': '4/mmm', '104': '4/mmm',
                        '105': '4/mmm', '106': '4/mmm', '107': '4/mmm', '108': '4/mmm', '109': '4/mmm', '110': '4/mmm',
                        '111': '4/mmm', '112': '4/mmm', '113': '4/mmm', '114': '4/mmm', '115': '4/mmm', '116': '4/mmm',
                        '117': '4/mmm', '118': '4/mmm', '119': '4/mmm', '120': '4/mmm', '121': '4/mmm', '122': '4/mmm',
                        '123': '4/mmm', '124': '4/mmm', '125': '4/mmm', '126': '4/mmm', '127': '4/mmm', '128': '4/mmm',
                        '129': '4/mmm', '130': '4/mmm', '131': '4/mmm', '132': '4/mmm', '133': '4/mmm', '134': '4/mmm',
                        '135': '4/mmm', '136': '4/mmm', '137': '4/mmm', '138': '4/mmm', '139': '4/mmm', '140': '4/mmm',
                        '141': '4/mmm', '142': '4/mmm', '143': '4/mmm', '144': '-3', '145': '-3', '146': '-3',
                        '147': '-3', '148': '-3', '149': '-3', '150': '-3m', '151': '-3m', '152': '-3m', '153': '-3m',
                        '154': '-3m', '155': '-3m', '156': '-3m', '157': '-3m', '158': '-3m', '159': '-3m',
                        '160': '-3m', '161': '-3m', '162': '-3m', '163': '-3m', '164': '-3m', '165': '-3m',
                        '166': '-3m', '167': '-3m', '168': '-3m', '169': '6/m', '170': '6/m', '171': '6/m',
                        '172': '6/m', '173': '6/m', '174': '6/m', '175': '6/m', '176': '6/m', '177': '6/m',
                        '178': '6/mmm', '179': '6/mmm', '180': '6/mmm', '181': '6/mmm', '182': '6/mmm', '183': '6/mmm',
                        '184': '6/mmm', '185': '6/mmm', '186': '6/mmm', '187': '6/mmm', '188': '6/mmm', '189': '6/mmm',
                        '190': '6/mmm', '191': '6/mmm', '192': '6/mmm', '193': '6/mmm', '194': '6/mmm', '195': '6/mmm',
                        '196': 'm-3', '197': 'm-3', '198': 'm-3', '199': 'm-3', '200': 'm-3', '201': 'm-3',
                        '202': 'm-3', '203': 'm-3', '204': 'm-3', '205': 'm-3', '206': 'm-3', '207': 'm-3',
                        '208': 'm-3m', '209': 'm-3m', '210': 'm-3m', '211': 'm-3m', '212': 'm-3m', '213': 'm-3m',
                        '214': 'm-3m', '215': 'm-3m', '216': 'm-3m', '217': 'm-3m', '218': 'm-3m', '219': 'm-3m',
                        '220': 'm-3m', '221': 'm-3m', '222': 'm-3m', '223': 'm-3m', '224': 'm-3m', '225': 'm-3m',
                        '226': 'm-3m', '227': 'm-3m', '228': 'm-3m', '229': 'm-3m', '230': 'm-3m'},
    }
    narrowed_down_notations = []                                    # Initialize list for tracking matching notations
    common_notation = True                                          # Assume notation is uniform until proved otherwise
    counter = 0                                                     # Initialize while loop counter
    total_phases = len(phases["point_group"])                       # Get number of phases
    if total_phases > 1:                                            # If more than one phase
        for phase in phases["point_group"]:                         # Check all header phases
            for notation in symmetry_notations:                     # Check all notation standards
                if phase not in laue_group_dict[notation]:          # If symmetry not found in current notation and
                    if notation in narrowed_down_notations:         # the notation that doesn't match is in the narrowed
                        narrowed_down_notations.remove(notation)    # down notations, eliminate as a possibility
                if phase in laue_group_dict[notation]:              # Is the current phase found in the current notation
                    if notation not in narrowed_down_notations:     # Is the matching notation in the narrowed_down list
                        narrowed_down_notations.append(notation)    # Add notation as possible notation for other phases
        while len(narrowed_down_notations) != 1 and counter > 10:   # Loop until only one notation possibility remaining
            counter += 1
            for phase in phases["point_group"]:
                for notation in narrowed_down_notations:            # Now go back through the possible assumed notations
                    if phase not in laue_group_dict[notation]:
                        if notation in narrowed_down_notations:
                            narrowed_down_notations.remove(notation)
                    if phase in laue_group_dict[notation]:
                        if notation not in narrowed_down_notations:
                            narrowed_down_notations.append(notation)
    else:
        for i, phase in enumerate(phases["point_group"]):                       # If only one phase
            for notation in symmetry_notations:
                if phase not in laue_group_dict[notation] and phase != '0':
                    warnings.warn(f"Input header symmetry unrecognized")        # User input symmetry is not recognized
    if len(narrowed_down_notations) != 1:                                       # Check to ensure uniform notation usage
        common_notation = False
    if common_notation is False:                                                # Warn user regarding non-uniformity
        warnings.warn(f"Input header symmetries are not identified as the same notation. "
                      f"Possible symmetry notations determined are {narrowed_down_notations}")
    else:
        for i, phase in enumerate(phases["point_group"]):                               # Loop through phases to convert
            phases["point_group"][i] = laue_group_dict[narrowed_down_notations][phase]  # Convert phases symmetries

    return phases["point_group"]


def file_writer(
    filename,
    xmap,
    index=None,
    image_quality_prop=None,
    confidence_index_prop=None,
    detector_signal_prop=None,
    pattern_fit_prop=None,
    extra_prop=None,
):
    """Write a :class:`~orix.crystal_map.crystal_map.CrystalMap` to an
    .ang file readable (at least) by MTEX and EDAX TSL OIM Analysis v7.

    The columns are phi1, Phi, phi2, x, y, image_quality,
    confidence_index, phase_id, detector_signal, and pattern_fit.

    Parameters in masked out or non-indexed points are set to

    * euler angles = 4 * pi

    * image quality = 0

    * confidence index = -1

    * phase ID = 0 if single phase or -1 if multi phase

    * pattern fit = 180

    * detector signal = 0

    * extra properties = 0

    Parameters
    ----------
    filename : str
        File name with an ".ang" file extension to write to.
    xmap : CrystalMap
        Crystal map to write to file.
    index : int, optional
        If the crystal map has more than one rotation/match and phase
        ID per point, this index can be set to write that "layer" of
        data to the file. For properties to be written as well, these
        must also have multiple values per point. To get the best match
        at every point, use None (default).
    image_quality_prop : str, optional
        Which map property to use as the image quality. If None
        (default), "iq" or "imagequality", if present, is used, or just
        zeros. If the property has more than one value per point and
        `index` is not given, only the first value is used.
    confidence_index_prop : str, optional
        Which map property to use as the confidence index. If None
        (default), "ci", "confidenceindex", "scores", or "correlation",
        if present, is used, or just zeros. If the property has more
        than one value per point and `index` is not given, only the
        first value is used.
    detector_signal_prop : str, optional
        Which map property to use as the detector signal. If None
        (default), "ds", or "detector_signal", if present, is used, or
        just zeros. If the property has more than one value per point
        and `index` is not given, only the first value is used.
    pattern_fit_prop : str, optional
        Which map property to use as the pattern fit. If None
        (default), "fit" or "patternfit", if present, is used, or just
        zeros. If the property has more than one value per point and
        `index` is not given, only the first value is used.
    extra_prop : str or list of str, optional
        One or multiple properties to add as extra columns in the .ang
        file, as a string or a list of strings. If None (default), no
        extra properties are added. If a property has more than one
        value per point and `index` is not given, only the first value
        is used.
    """
    if xmap.ndim > 2:
        raise ValueError("Writing a 3D dataset to an .ang file is not supported")
    header = _get_header_from_phases(xmap)

    # Number of decimals to round to
    decimals = 5

    # Get map data, accounting for potentially masked out values
    nrows, ncols, dy, dx = _get_nrows_ncols_step_sizes(xmap)
    map_size = nrows * ncols
    # Rotations
    if index is not None:
        eulers = xmap.get_map_data(
            xmap.rotations[:, index].to_euler(), decimals=decimals, fill_value=0
        )
    else:
        eulers = xmap.get_map_data("rotations", decimals=decimals, fill_value=0)
    eulers = eulers.reshape((map_size, 3))
    indexed_points = xmap.get_map_data(xmap.is_indexed, fill_value=False).reshape(
        map_size
    )
    eulers[~indexed_points] = 4 * np.pi

    # Coordinate arrays
    d, _ = create_coordinate_arrays((nrows, ncols), (dy, dx))
    x = d["x"]
    y = d["y"]
    x_width = _get_column_width(np.max(x))
    y_width = _get_column_width(np.max(y))

    # Properties
    desired_prop_names = [
        image_quality_prop,
        confidence_index_prop,
        detector_signal_prop,
        pattern_fit_prop,
    ]
    if extra_prop is not None:
        if isinstance(extra_prop, str):
            extra_prop = [
                extra_prop,
            ]
        desired_prop_names += extra_prop
    prop_arrays = _get_prop_arrays(
        xmap=xmap,
        prop_names=list(xmap.prop.keys()),
        desired_prop_names=desired_prop_names,
        map_size=map_size,
        index=index,
        decimals=decimals,
    )
    # Set values for non-indexed points
    prop_arrays[~indexed_points, 0::2] = 0  # IQ, detector signal
    prop_arrays[~indexed_points, 1] = -1  # CI
    prop_arrays[~indexed_points, 3] = 180  # Pattern fit
    prop_arrays[~indexed_points, 4:] = 0
    # Get property column widths
    prop_widths = [
        _get_column_width(np.max(prop_arrays[:, i]))
        for i in range(prop_arrays.shape[1])
    ]

    # Phase ID
    original_phase_ids = xmap.get_map_data("phase_id").reshape(map_size)
    pl = xmap.phases.deepcopy()
    if -1 in pl.ids:
        del pl[-1]
    if pl.size > 1:
        new_phase_ids = np.zeros(map_size, dtype=int)
    else:
        new_phase_ids = -np.ones(map_size, dtype=int)
    # Phase IDs are reversed because EDAX TSL OIM Analysis v7.2.0
    # assumes a reversed phase order in the header
    for i, (phase_id, phase) in reversed(list(enumerate(pl))):
        new_phase_ids[original_phase_ids == phase_id] = i + 1

    # Extend header with column names
    header += (
        "\n"
        "Column names: phi1, Phi, phi2, x, y, image_quality, confidence_index, "
        "phase_id, detector_signal, pattern_fit"
    )

    # Get data formats
    fmt = (
        f"%8.5f  %8.5f  %8.5f  %{x_width}.5f  %{y_width}.5f  %{prop_widths[0]}.5f  "
        f"%{prop_widths[1]}.5f  %2i  %{prop_widths[2]}.5f  %{prop_widths[3]}.5f"
    )
    if extra_prop is not None:  # Then it must be a list!
        for extra_width, extra_prop_name in zip(prop_widths[4:], extra_prop):
            fmt += f"  %{extra_width}.5f"
            header += f", {extra_prop_name}"
    header += "\n"

    # Finally, write everything to file
    np.savetxt(
        fname=filename,
        X=np.column_stack(
            [eulers, x, y, prop_arrays[:, :2], new_phase_ids, prop_arrays[:, 2:]]
        ),
        fmt=fmt,
        header=header,
    )


def _get_header_from_phases(xmap):
    """Return a string with the .ang file header from the crystal map
    metadata.

    Notes
    -----
    Phase IDs are reversed because EDAX TSL OIM Analysis v7.2.0 assumes
    a reversed phase order in the header
    """
    nrows, ncols, _, _ = _get_nrows_ncols_step_sizes(xmap)
    # Initialize header with microscope info
    header = (
        "TEM_PIXperUM           1.000000\n"
        "x-star                 0.000000\n"
        "y-star                 0.000000\n"
        "z-star                 0.000000\n"
        "WorkingDistance        0.000000\n"
        "\n"
    )
    # Get phase list, removing the non indexed phase if present
    pl = xmap.phases.deepcopy()
    if -1 in pl.ids:
        del pl[-1]

    # Extend header with phase info
    # Phase IDs are reversed because EDAX TSL OIM Analysis v7.2.0
    # assumes a reversed phase order in the header
    for i, (_, phase) in reversed(list(enumerate(pl))):
        lattice_constants = phase.structure.lattice.abcABG()
        lattice_constants = " ".join([f"{float(val):.3f}" for val in lattice_constants])
        phase_id = i + 1
        phase_name = phase.name
        if phase_name == "":
            phase_name = f"phase{phase_id}"
        if phase.point_group is None:
            point_group_name = "1"
        else:
            proper_point_group = phase.point_group.proper_subgroup
            point_group_name = proper_point_group.name
            for key, alias in point_group_aliases.items():
                if point_group_name == key:
                    point_group_name = alias[0]
                    break
        header += (
            f"Phase {phase_id}\n"
            f"MaterialName    {phase_name}\n"
            f"Formula    {phase_name}\n"
            f"Info\n"
            f"Symmetry    {point_group_name}\n"
            f"LatticeConstants    {lattice_constants}\n"
            "NumberFamilies    0\n"
        )
    # Extend header with map info
    header += (
        "GRID: SqrGrid\n"
        f"XSTEP: {float(xmap.dx):.6f}\n"
        f"YSTEP: {float(xmap.dy):.6f}\n"
        f"NCOLS_ODD: {ncols}\n"
        f"NCOLS_EVEN: {ncols}\n"
        f"NROWS: {nrows}\n"
        "\n"
        f"OPERATOR: orix_v{__version__}\n"
        "\n"
        "SAMPLEID:\n"
        "\n"
        "SCANID:\n"
    )
    return header


def _get_nrows_ncols_step_sizes(xmap):
    """Get crystal map shape and step sizes."""
    nrows, ncols = (1, 1)
    dy, dx = xmap.dy, xmap.dx
    if xmap.ndim == 1:
        ncols = xmap.shape[0]
        dy = 1
    else:  # xmap.ndim == 2:
        nrows, ncols = xmap.shape
    return nrows, ncols, dy, dx


def _get_column_width(max_value, decimals=5):
    """Get width of column to pass to :func:`numpy.savetxt`, accounting
    for the decimal point and a sign +/-.
    """
    return len(str(int(max_value // 1))) + decimals + 2


def _get_prop_arrays(xmap, prop_names, desired_prop_names, map_size, index, decimals=5):
    """Return a 2D array (n_points, n_properties) with desired property
    values in, or just zeros.

    This function tries to get as many properties as possible from the
    crystal map properties.

    Parameters
    ----------
    xmap : CrystalMap
    prop_names : list of str
    desired_prop_names : list of str
    map_size : int
    index : int or None
    decimals : int, optional

    Returns
    -------
    numpy.ndarray
    """
    # "Image_quality" -> "imagequality" etc.
    prop_names_lower = [k.lower().replace("_", "") for k in prop_names]
    prop_names_lower_arr = np.array(prop_names_lower)
    # Potential extra names added so that lists are of the same length
    # in the loop
    all_expected_prop_names = [
        ["iq", "imagequality"],
        ["ci", "confidenceindex", "scores", "correlation"],
        ["ss", "semsignal", "detectorsignal"],
        ["fit", "patternfit"],
    ] + desired_prop_names[4:]
    n_desired_props = len(desired_prop_names)
    prop_arrays = np.zeros((map_size, n_desired_props), dtype=np.float32)
    for i, (name, names) in enumerate(zip(desired_prop_names, all_expected_prop_names)):
        prop = _get_prop_array(
            xmap=xmap,
            prop_name=name,
            expected_prop_names=names,
            prop_names=prop_names,
            prop_names_lower_arr=prop_names_lower_arr,
            decimals=decimals,
            index=index,
            fill_value=0,
        )
        if prop is not None:
            prop_arrays[:, i] = prop.reshape(map_size)
    return prop_arrays


def _get_prop_array(
    xmap,
    prop_name,
    expected_prop_names,
    prop_names,
    prop_names_lower_arr,
    index,
    decimals=5,
    fill_value=0,
):
    """Return a 1D array (n_points,) with the desired property values or
    None if the property cannot be read.

    Reasons for why the property cannot be read:
    * Property name isn't among the crystal map properties
    * Property has only one value per point, but `index` is not None

    Parameters
    ----------
    xmap : CrystalMap
    prop_name : str
    expected_prop_names : list of str
    prop_names : list of str
    prop_names : list of str
    index : int or None
    decimals : int, optional
    fill_value : int, float, or bool, optional

    Returns
    -------
    numpy.ndarray or None
    """
    kwargs = dict(decimals=decimals, fill_value=fill_value)
    if len(prop_names_lower_arr) == 0 and prop_name is None:
        return
    else:
        if prop_name is None:
            # Search for a suitable property
            for k in expected_prop_names:
                is_equal = k == prop_names_lower_arr
                if is_equal.any():
                    prop_name = prop_names[np.argmax(is_equal)]
                    break
            else:  # If no suitable property was found
                return
        # There is a property
        if len(xmap.prop[prop_name].shape) == 1:
            # Return the single array even if `index` is given
            return xmap.get_map_data(prop_name, **kwargs)
        else:
            if index is None:
                index = 0
            return xmap.get_map_data(xmap.prop[prop_name][:, index], **kwargs)