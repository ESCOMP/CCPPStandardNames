#!/usr/bin/env python

"""
Remove duplicates from a metadata standard-name XML library file.
"""

import argparse
import sys
import os.path
import xml.etree.ElementTree as ET
from xml_tools import read_xml_file
import copy

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("standard_name_file",
                        metavar='<standard names filename>',
                        type=str, help="XML file with standard name library")
    pargs = parser.parse_args(args)
    return pargs

###############################################################################
def main_func():
###############################################################################
    """Parse the standard names database file and notify of duplicates.
    """
    # Parse command line arguments
    args = parse_command_line(sys.argv[1:], __doc__)
    stdname_file = os.path.abspath(args.standard_name_file)    
    _, root = read_xml_file(stdname_file)
    
    #get list of all standard names
    all_std_names = []
    for name in root.findall('./section/standard_name'):
        all_std_names.append(name.attrib['name'])
    
    #get list of all unique and duplicate standard names, in source order
    seen = set()
    uniq_std_names = []
    dup_std_names = []
    for x in all_std_names:
        if x not in seen:
            uniq_std_names.append(x)
            seen.add(x)
        else:
            dup_std_names.append(x)
    
    #delete all duplicate elements after the first
    for dup in dup_std_names:
        rm_elements = root.findall('./section/standard_name[@name="%s"]'%dup)[1:]
        print("{0}, ({1} duplicate(s))".format(dup, len(rm_elements)))
        
        #TODO:
        #rm_parents = root.findall('./section/standard_name[@name="%s"]...'%dup)[1:]
        #for i in range(len(rm_elements)):
        #    print(ET.tostring(rm_elements[i]))
        #    rm_parents[i].remove(rm_elements[i])
    
    #TODO: overwrite XML file with duplicates removed
    
###############################################################################
if __name__ == "__main__":
    main_func()