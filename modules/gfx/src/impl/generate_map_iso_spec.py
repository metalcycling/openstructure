#------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------
#
# Author: Ansgar Philippsen
#
# automatically generate iso-contouring marching cube specs
#

import math

# cube definition

# corners:
#     5        6
#      -------
#     /.     /|
#    / .    / |
#  1 ------ 2 |
#   |  ....|..|
#   | .4   | / 7
#   |.     |/
#    ------
#  0        3

# vertices:              
#      ---9---
#     5.     6|
#    / 8    / |
#    ---1--   10
#   |  ..11|..|
#   0 4    2 7  
#   |.     |/
#    --3---

# faces     2
#          .       
#      -------
#     /      /|
#    /   1  / |
#    ------  .|.. 5
#4..|      |  |
#   |  0   | / 
#   |      |/.
#    ------   3

corner_coord=[[0,0,0],[0,1,0],[1,1,0],[1,0,0],
               [0,0,1],[0,1,1],[1,1,1],[1,0,1]]

edge_corner=[[0,1],[1,2],[2,3],[3,0],
                [0,4],[1,5],[2,6],[3,7],
                [4,5],[5,6],[6,7],[4,7]]

face_edges=[[0,1,2,3],[1,5,9,6],[9,8,11,10],
            [11,7,3,4],[4,8,5,0],[7,2,6,10]]

header_code="""
// automatically generated by generate_map_iso_spec.py
// do not edit directly

#include "map_iso_gen.hh"

namespace ost { namespace gfx { namespace map_iso {

template<int N>
void AddLinesAndFaces(IndexedVertexArray& va, // vertex array
		      unsigned int vertex_id[12] // this list of vertex ids
		      );
"""

isocube_header_code="""
template<> 
void AddLinesAndFaces<%d>(IndexedVertexArray& va,unsigned int vertex_id[12] )
{
"""

isocube_addline_code="""  va.AddLine(vertex_id[%d],vertex_id[%d]);
"""

isocube_addface_code="""  va.AddTriN(vertex_id[%d],vertex_id[%d],vertex_id[%d]);
"""

isocube_debug_code="""
#ifdef MAP_ISO_DEBUG
  if(vertex_id[%d]==0) std::cerr << \"ISO_DEBUG: \" << %d << \" \" << %d << std::endl;
  if(vertex_id[%d]==0) std::cerr << \"ISO_DEBUG: \" << %d << \" \" << %d << std::endl;
  if(vertex_id[%d]==0) std::cerr << \"ISO_DEBUG: \" << %d << \" \" << %d << std::endl;
#endif
"""

isocube_footer_code="}\n"

flist_header_code="""
void GenerateLFList(AddLFList& lfl)
{
  lfl.clear();
"""

flist_code="  lfl.push_back(AddLinesAndFaces<%d>);\n"

flist_footer_code="""
}
"""

footer_code="""
}}} // ns
"""

def vadd(v1,v2):
    return [v1[0]+v2[0],
            v1[1]+v2[1],
            v1[2]+v2[2]]

def vsub(v1,v2):
    return [v1[0]-v2[0],
            v1[1]-v2[1],
            v1[2]-v2[2]]

def vmul(s,v1):
    return [v1[0]*s,
            v1[1]*s,
            v1[2]*s]

def vdot(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def vcross(v1,v2):
    return [v1[1]*v2[2]-v2[1]*v1[2],
	    v1[2]*v2[0]-v2[2]*v1[0],
	    v1[0]*v2[1]-v2[0]*v1[1]]

def vnorm(v1):
    len2=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]
    return vmul(1.0/math.sqrt(len2),v1)

def work(spec_file,pattern):

    corner_bit=[0,0,0,0,0,0,0,0]
    for i in range(8):
        if pattern & (1<<i):
            corner_bit[i]=1

    edge_bit=[0,0,0,0,0,0,0,0,0,0,0,0]
    edge_list=[]
    for i in range(12):
        if corner_bit[edge_corner[i][0]] != corner_bit[edge_corner[i][1]]:
            edge_bit[i]=1
            edge_list.append(i)

    spec_file.write(isocube_header_code%pattern)

    # first the lines
    edge_pairs=[]
    for e1 in edge_list:
        for e2 in edge_list:
            if e1>e2:
                for f in face_edges:
                    if e1 in f and e2 in f:
                        # line across a face
                        edge_pairs.append([e1,e2])
                        spec_file.write(isocube_addline_code%(e1,e2))

    # find line-loop(s)
    while len(edge_pairs)>0:
        edge_loop=[]
        edge_loop.extend(edge_pairs[0])
        edge_pairs.remove(edge_pairs[0])
        closed=False
        while not closed:
            for ep in edge_pairs:
                if ep[0]==edge_loop[len(edge_loop)-1]:
                    if ep[1]==edge_loop[0]:
                        closed=True
                    else:
                        edge_loop.append(ep[1])
                    edge_pairs.remove(ep)
                    break
                elif ep[1]==edge_loop[len(edge_loop)-1]:
                    if ep[0]==edge_loop[0]:
                        closed=True
                    else:
                        edge_loop.append(ep[0])
                    edge_pairs.remove(ep)
                    break

        # find correct orientation 
        # assemble dir0 from middle of edge and the higher value corner
        edge_id0 = edge_loop[0]
        v00 = edge_corner[edge_id0][0]
        v01 = edge_corner[edge_id0][1]
        e0_pos = vmul(0.5,vadd(corner_coord[v00],corner_coord[v01]))
        e1_pos=[0,0,0]
        if corner_bit[v00]==1:
            e1_pos=corner_coord[v00]
        else:
            e1_pos=corner_coord[v01]
        dir0 = vnorm(vsub(e1_pos,e0_pos))
        # assemble dir1 from the cross product of the lines of the
        # edge_loop
        edge_id1 = edge_loop[1]
        v10 = edge_corner[edge_id1][0]
        v11 = edge_corner[edge_id1][1]
        edge_id2 = edge_loop[2]
        v20 = edge_corner[edge_id2][0]
        v21 = edge_corner[edge_id2][1]
        e1_pos = vmul(0.5,vadd(corner_coord[v10],corner_coord[v11]))
        e2_pos = vmul(0.5,vadd(corner_coord[v20],corner_coord[v21]))
        e10 = vsub(e0_pos,e1_pos)
        e12 = vsub(e2_pos,e1_pos)
        dir1 = vnorm(vcross(e10,e12))
        # this is the actual orientation check
        if vdot(dir0,dir1)<0.0:
            edge_loop.reverse()
        # dump code to file
        for i in range(len(edge_loop)-2):
            spec_file.write(isocube_addface_code%(edge_loop[0],edge_loop[i+2],edge_loop[i+1]))
            spec_file.write(isocube_debug_code%(edge_loop[0],edge_loop[0],pattern,
                                                edge_loop[i+2],edge_loop[i+2],pattern,
                                                edge_loop[i+1],edge_loop[i+1],pattern))
                                                
            
    spec_file.write(isocube_footer_code)
    

spec_file=file("map_iso_spec.hh","w")

spec_file.write(header_code)
    
for i in range(256):
    work(spec_file,i)

spec_file.write(flist_header_code)
for i in range(256):
    spec_file.write(flist_code%i)
spec_file.write(flist_footer_code)

spec_file.write(footer_code)