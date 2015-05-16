#!/usr/bin/python3

import os;
import csv;
import argparse;

parser = argparse.ArgumentParser();
parser.add_argument("folder", help="Folder to start from")
args = parser.parse_args();
if not args.folder:
	os.exit(1);

#root = "./png_noise4_n0.05";
root = args.folder;
reader = csv.reader(open("%s/stats.csv"%root, mode='r'))
mydict = dict((rows[0],rows[1]) for rows in reader)
entries = {};
for item in os.listdir(root):
	if ( item.find("png") >= 0 ):
		keyStr = item[ item.find("pw")+2 : item.find(".png") ];
		#print( keyStr );
		entries[ float(keyStr) ] = [ keyStr, item, mydict[keyStr] ];

print( "<html><body>" );
for key in sorted(entries):
	print( "<div style=\"float:left\">");
	print( "<img src=\"%s\" width=\"320\">" % (entries[key][1]) );
	print( "<p align=\"center\">%s,#%s</p>" % (entries[key][0],entries[key][2]) );
	print( "</div>");

print("<div style=\"clear:both\" />");

print( "<table><tr>");
for key in sorted(entries):
	print("<th>%s</th>" % (entries[key][0]) );
print("</tr><tr>");
for key in sorted(entries):
	print("<td>%s</td>" % (entries[key][2]) );
print("</tr>");
print("</body></html>");