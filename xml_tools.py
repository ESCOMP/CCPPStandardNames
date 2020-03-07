#!/usr/bin/env python

"""
Parse and / or validate an XML file and return the captured variables.
"""

# Python library imports
from __future__ import print_function
import os
import os.path
import subprocess
import sys
import logging
from distutils.spawn import find_executable
import xml.etree.ElementTree as ET
try:
    _XMLLINT = find_executable('xmllint')
except ImportError:
    _XMLLINT = None
# end try

# Find python version
PY3 = sys.version_info[0] > 2
PYSUBVER = sys.version_info[1]
_LOGGER = None

###############################################################################
def call_command(commands, logger, silent=False):
###############################################################################
    """
    Try a command line and return the output on success (None on failure)
    >>> call_command(['ls', 'really__improbable_fffilename.foo'], _LOGGER) #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    RuntimeError: Execution of 'ls really__improbable_fffilename.foo' failed:
    [Errno 2] No such file or directory
    >>> call_command(['ls', 'really__improbable_fffilename.foo'], _LOGGER, silent=True)
    False
    >>> call_command(['ls'], _LOGGER)
    True
    """
    result = False
    outstr = ''
    if logger is None:
        silent = True
    # end if
    try:
        if PY3:
            if PYSUBVER > 6:
                cproc = subprocess.run(commands, check=True,
                                       capture_output=True,
                                       stderr=subprocess.STDOUT)
                if not silent:
                    logger.debug(cproc.stdout)
                # end if
                result = cproc.returncode == 0
            elif PYSUBVER >= 5:
                cproc = subprocess.run(commands, check=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
                if not silent:
                    logger.debug(cproc.stdout)
                # end if
                result = cproc.returncode == 0
            else:
                raise ValueError("Python 3 must be at least version 3.5")
            # end if
        else:
            pproc = subprocess.Popen(commands, stdin=None,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            output, _ = pproc.communicate()
            if not silent:
                logger.debug(output)
            # end if
            result = pproc.returncode == 0
        # end if
    except (OSError, RuntimeError, subprocess.CalledProcessError) as err:
        if silent:
            result = False
        else:
            cmd = ' '.join(commands)
            emsg = "Execution of '{}' failed with code:\n"
            outstr = emsg.format(cmd, err.returncode)
            outstr += "{}".format(err.output)
            raise RuntimeError(outstr)
        # end if
    # end of try
    return result

###############################################################################
def find_schema_version(root):
###############################################################################
    """
    Find the version of the host registry file represented by root
    >>> find_schema_version(ET.fromstring('<model name="CAM" version="1.0"></model>'))
    [1, 0]
    >>> find_schema_version(ET.fromstring('<model name="CAM" version="1.a"></model>')) #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: Illegal version string, '1.a'
    Format must be <integer>.<integer>
    >>> find_schema_version(ET.fromstring('<model name="CAM" version="0.0"></model>')) #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: Illegal version string, '0.0'
    Major version must be at least 1
    >>> find_schema_version(ET.fromstring('<model name="CAM" version="0.-1"></model>')) #doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: Illegal version string, '0.0'
    Minor version must be at least 0
    """
    verbits = None
    if 'version' not in root.attrib:
        raise ValueError("version attribute required")
    # end if
    version = root.attrib['version']
    versplit = version.split('.')
    try:
        if len(versplit) != 2:
            raise ValueError('oops')
        # end if (no else needed)
        try:
            verbits = [int(x) for x in versplit]
        except ValueError as verr:
            raise ValueError(verr)
        # end try
        if verbits[0] < 1:
            raise ValueError('Major version must be at least 1')
        # end if
        if verbits[1] < 0:
            raise ValueError('Minor version must be non-negative')
        # end if
    except ValueError as verr:
        errstr = """Illegal version string, '{}'
        Format must be <integer>.<integer>"""
        ve_str = str(verr)
        if ve_str:
            errstr = ve_str + '\n' + errstr
        # end if
        raise ValueError(errstr.format(version))
    # end try
    return verbits

###############################################################################
def find_schema_file(schema_root, version, schema_path=None):
###############################################################################
    """Find and return the schema file based on <schema_root> and <version>
    or return None.
    If <schema_path> is present, use that as the directory to find the
    appropriate schema file. Otherwise, just look in the current directory."""

    verstring = '_'.join([str(x) for x in version])
    schema_filename = "{}_v{}.xsd".format(schema_root, verstring)
    if schema_path:
        schema_file = os.path.join(schema_path, schema_filename)
    else:
        schema_file = schema_filename
    # end if
    if os.path.exists(schema_file):
        return schema_file
    # end if
    return None

###############################################################################
def validate_xml_file(filename, schema_root, version, logger,
                      schema_path=None, error_on_noxmllint=False):
###############################################################################
    """
    Find the appropriate schema and validate the XML file, <filename>,
    against it using xmllint
    """
    # Check the filename
    if not os.path.isfile(filename):
        raise ValueError("validate_xml_file: Filename, '{}', does not exist".format(filename))
    # end if
    if not os.access(filename, os.R_OK):
        raise ValueError("validate_xml_file: Cannot open '{}'".format(filename))
    # end if
    if not schema_path:
        # Find the schema, based on the model version
        thispath = os.path.abspath(__file__)
        pdir = os.path.dirname(os.path.dirname(os.path.dirname(thispath)))
        schema_path = os.path.join(pdir, 'schema')
    # end if
    schema_file = find_schema_file(schema_root, version, schema_path)
    if not (schema_file and os.path.isfile(schema_file)):
        verstring = '.'.join([str(x) for x in version])
        emsg = """validate_xml_file: Cannot find schema for version {},
        {} does not exist"""
        raise ValueError(emsg.format(verstring, schema_file))
    # end if
    if not os.access(schema_file, os.R_OK):
        emsg = "validate_xml_file: Cannot open schema, '{}'"
        raise ValueError(emsg.format(schema_file))
    # end if
    if _XMLLINT is not None:
        if logger is not None:
            lmsg = "Checking file {} against schema {}"
            logger.debug(lmsg.format(filename, schema_file))
        # end if
        cmd = [_XMLLINT, '--noout', '--schema', schema_file, filename]
        result = call_command(cmd, logger)
        return result
    # end if
    lmsg = "xmllint not found, could not validate file {}"
    if error_on_noxmllint:
        raise ValueError("validate_xml_file: " + lmsg.format(filename))
    # end if
    if logger is not None:
        logger.warning(lmsg.format(filename))
    # end if
    return True # We could not check but still need to proceed

###############################################################################
def read_xml_file(filename, logger=None):
###############################################################################
    """Read the XML file, <filename>, and return its tree and root"""
    if os.path.isfile(filename) and os.access(filename, os.R_OK):
        if PY3:
            file_open = (lambda x: open(x, 'r', encoding='utf-8'))
        else:
            file_open = (lambda x: open(x, 'r'))
        # end if
        with file_open(filename) as file_:
            try:
                tree = ET.parse(file_)
                root = tree.getroot()
            except ET.ParseError as perr:
                emsg = "read_xml_file: Cannot read {}, {}"
                raise ValueError(emsg.format(filename, perr))
    elif not os.access(filename, os.R_OK):
        raise ValueError("read_xml_file: Cannot open '{}'".format(filename))
    else:
        emsg = "read_xml_file: Filename, '{}', does not exist"
        raise ValueError(emsg.format(filename))
    # end if
    if logger:
        logger.debug("Read XML file, '{}'".format(filename))
    # end if
    return tree, root

###############################################################################

if __name__ == "__main__":
    _LOGGER = logging.getLogger('xml_tools')
    for handler in list(_LOGGER.handlers):
        _LOGGER.removeHandler(handler)
    # end for
    _LOGGER.addHandler(logging.NullHandler())
    try:
        # First, run doctest
        import doctest
        doctest.testmod()
    except ValueError as cerr:
        print("{}".format(cerr))
# no else:
