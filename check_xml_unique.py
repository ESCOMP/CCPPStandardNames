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
    parser.add_argument("--overwrite", action='store_true',
                        help="flag to remove duplicates and overwrite the file")
    
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
    
    if args.overwrite:
        #delete all duplicate elements after the first
        if len(dup_std_names)>0:
            print('The following duplicate standard names were found:')
            for dup in dup_std_names:
                rm_elements = root.findall('./section/standard_name[@name="%s"]'%dup)[1:]
                print("{0}, ({1} duplicate(s))".format(dup, len(rm_elements)))
            print('Removing duplicates and overwriting {}'.format(stdname_file))
            for dup in dup_std_names:
                rm_parents = root.findall('./section/standard_name[@name="%s"]...'%dup)[1:]
                for par in rm_parents:
                    rm_ele = par.findall('./standard_name[@name="%s"]'%dup)
                    for ele in rm_ele:
                        par.remove(ele)
            
            _.write(stdname_file, "utf-8")
        else:
            print('No duplicate standard names were found.')
    else:
        #write out duplicate standard names
        if len(dup_std_names)>0:
            print('The following duplicate standard names were found:')
            for dup in dup_std_names:
                rm_elements = root.findall('./section/standard_name[@name="%s"]'%dup)[1:]
                print("{0}, ({1} duplicate(s))".format(dup, len(rm_elements)))
        else:
            print('No duplicate standard names were found.')
    
###############################################################################
if __name__ == "__main__":
    main_func()