#!/usr/bin/env python3

"""
Remove duplicates from a metadata standard-name XML library file.
"""

import argparse
import sys
import os.path
import copy
import re

################################################
#Add CCPP framework (lib) modules to python path
################################################

_CURR_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(_CURR_DIR, "lib"))

#######################################
#Import needed framework python modules
#######################################

from xml_tools import find_schema_file, find_schema_version, validate_xml_file, read_xml_file

def main():
    """Parse the standard names database file and output a dictionary
    where the keys are any standard names in violation of character rules,
    and the values are lists of the specific rules violated
    """
    #Parse arguments
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-s","--standard_name_file",
                        metavar='<standard names filename>',required=True,
                        type=str, help="XML file with standard name library")
    args = parser.parse_args()

    stdname_file = os.path.abspath(args.standard_name_file)
    tree, root = read_xml_file(stdname_file)

    # Validate the XML file
    version = find_schema_version(root)
    schema_name = os.path.basename(stdname_file)[0:-4]
    schema_root = os.path.dirname(stdname_file)
    schema_path = os.path.join(schema_root,schema_name)
    schema_file = find_schema_file(schema_path, version)
    if schema_file:
        try:
            validate_xml_file(stdname_file, schema_name, version, None,
                            schema_path=schema_root, error_on_noxmllint=True)
        except ValueError:
            raise ValueError("Invalid standard names file, {}".format(stdname_file))
    else:
        raise ValueError(
            'Cannot find schema file, {}, for version {}'.format(schema_name, version)
        )

    #Parse list of standard names and see if any names violate one or more rules
    violators = {}
    for name in root.findall('./section/standard_name'):
        sname = name.attrib['name']
        violations = []
        if re.sub('[a-z]', '', sname[0]):
            violations.append('First character is not a lowercase letter')
        testchars = re.sub('[a-z0-9_]', '', sname)
        if testchars:
            violations.append(f'Invalid characters are present: {testchars}')

        # If any violations were detected, add an entry to "violators" dictionary
        if violations:
            violators[sname] = violations

    if violators:
        raise Exception(f"Violating standard names found:\n{violators}")
    else:
        print(f'Success! All standard names in {args.standard_name_file} follow the rules.')

if __name__ == "__main__":
    main()
