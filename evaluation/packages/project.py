"""@package docstring
This module defines an interface to read and use InputGen projects

See C++ InputGen project for more details on Projects
"""

import xml.etree.ElementTree as ET

class PyProject:
    """Unique class of the module
    """

    def __init__(self, path):
        """Default constructor, taking as input a xml file."""
        print 'Loading project ' + path.name

        tree = ET.parse(path)
        root = tree.getroot()
        for child in root:
            print child.tag, child.attrib

