#!/usr/bin/env python

import math
from pyacs.lib.euclid import *

def distance(A,B):
    return(math.sqrt((A[0]-B[0])**2+(A[1]-B[1])**2+(A[2]-B[2])**2))

def icosahedron():
    """Constructs an icosahedron on the unit-sphere"""

    tau= (1 + math.sqrt(5)) / 2
    # normalized
    radius=math.sqrt(1+tau**2)

    verts = [ \
        (0,  1,  tau),
        (0,  1, -tau),
        (0, -1,  tau),
        (0, -1, -tau),

        ( 1,  tau, 0),
        ( 1, -tau, 0),
        (-1,  tau, 0),
        (-1, -tau, 0),
        
        ( tau, 0,  1),
        ( tau, 0, -1),
        (-tau, 0,  1),
        (-tau, 0, -1)]

    

    
    faces=[\
        (0, 2, 8),
        (0, 2, 10),
        (0, 4, 6),
        (0, 4, 8),
        (0, 6, 10),
        (1, 3, 9),
        (1, 3, 11),
        (1, 4, 6),
        (1, 4, 9),
        (1, 6, 11),
        (2, 5, 7),
        (2, 5, 8),
        (2, 7, 10),
        (3, 5, 7),
        (3, 5, 9),
        (3, 7, 11),
        (4, 8, 9),
        (5, 8, 9),
        (6, 10, 11),
        (7, 10, 11)]

    # normalize
    verts = [ \
        (0,  1/radius,  tau/radius),
        (0,  1/radius, -tau/radius),
        (0, -1/radius,  tau/radius),
        (0, -1/radius, -tau/radius),

        ( 1/radius,  tau/radius, 0),
        ( 1/radius, -tau/radius, 0),
        (-1/radius,  tau/radius, 0),
        (-1/radius, -tau/radius, 0),
        
        ( tau/radius, 0,  1/radius),
        ( tau/radius, 0, -1/radius),
        (-tau/radius, 0,  1/radius),
        (-tau/radius, 0, -1/radius)]
    return(verts,faces)

def subdivide(verts, faces):
    """Subdivide each triangle into four triangles, pushing verts to the unit sphere"""

    from progress.bar import Bar

    triangles = len(faces)
    bar = Bar('', max=triangles, suffix='%(percent).1f%% - %(elapsed)ds')
    for faceIndex in range(triangles):


        # Create three new verts at the midpoints of each edge:
        face = faces[faceIndex]
        a,b,c = (Vector3(*verts[vertIndex]) for vertIndex in face)
        verts.append((a + b).normalized()[:])
        verts.append((b + c).normalized()[:])
        verts.append((a + c).normalized()[:])

        # Split the current triangle into four smaller triangles:
        i = len(verts) - 3
        j, k = i+1, i+2
        faces.append((i, j, k))
        faces.append((face[0], i, k))
        faces.append((i, face[1], j))
        faces[faceIndex] = (k, j, face[2])

        bar.next()

    bar.finish()

    return verts, faces


def mesh_global(num_subdivisions=6):

    Rt=6371.0E3
    (verts,faces)=icosahedron()
    print("-- Number of vertices and faces of initial icosahedron ",len(verts),len(faces))
    print("-- Number of divisions ",num_subdivisions)
    print("-- Now doing subdivision...")
    for x in range(num_subdivisions):
        print("   - Division iteration: ",x+1,"/",num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ",len(verts),len(faces))
        A=verts[faces[0][0]]
        B=verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A,B)*Rt/1000.)))
    return (verts,faces)
    
    
    
    
def mesh_regional(num_subdivisions=6, bounds=None):
    """
    Makes an equilateral triangles mesh over a sphere of unit radius using successive subdivisions
    of an initial icosahedron 
    returns verts and faces
    verts: list of vertices ; a vertice is a list of [x,y,z] in geocentric cartesian coordinates over the sphere of unit radius
    faces: list of faces ; a face is a list of vertices index
    faces[i][j] gives the [X, Y, Z] coordinates of vertex j of face i
    """
    (lon_min,lon_max,lat_min,lat_max)=list(map(float,bounds.split('/')[1:]))
    (rlon_min,rlon_max,rlat_min,rlat_max)=list(map(math.radians,(lon_min,lon_max,lat_min,lat_max)))

    ###########################################################################
    # import
    ###########################################################################
 #   import Polygon
    from pyacs.lib import coordinates
    from shapely.geometry import Polygon as shp_polygon

 #   rectangle=Polygon.Polygon(((rlon_min,rlat_min), (rlon_max,rlat_min),(rlon_max,rlat_max),(rlon_min,rlat_max)))
    shp_rectangle =  shp_polygon( [(rlon_min,rlat_min), (rlon_max,rlat_min),(rlon_max,rlat_max),(rlon_min,rlat_max)] )

    Rt=6371.0E3
    (verts,faces)=icosahedron()
    print("-- Number of vertices and faces of initial icosahedron ",len(verts),len(faces))
    print("-- Number of divisions ",num_subdivisions)
    print("-- Now doing subdivision...")
    ndivision=6
    if num_subdivisions<6:ndivision=num_subdivisions
    for x in range(ndivision):
        print("   - Division iteration: ",x+1,"/",num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ",len(verts),len(faces))
        A=verts[faces[0][0]]
        B=verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A,B)*Rt/1000.)))
    
    print("-- Now doing selection on bounds ",bounds)
    lfaces=[]
    for face in faces:
        #print "#face: ",face
        lgeo=[]
        for j in [0,1,2]:
            (x,y,z)=verts[face[j]]
            (rlon,rlat,rh)=coordinates.xyz2geospheric(x,y,z)
            lgeo.append((rlon,rlat))
        
#        triangle=Polygon.Polygon((lgeo))
        shp_triangle = shp_polygon( lgeo )
#        intersection = rectangle & triangle
        shp_intersection = shp_rectangle.intersects( shp_triangle )
#        print('intersection ' , intersection)
#        print('shp_intersection ' , shp_intersection)
#        if not intersection:continue#print "No intersection"
        if not shp_intersection:continue#print "No intersection"
        else:
            #print "Intersection"    
            lfaces.append(face)
    
    print("-- Keeping ",len(lfaces)," faces")
    faces=lfaces
    
    print("-- Now dealing with finer subdivisions")
    for x in range(ndivision,num_subdivisions):
        print("   - Division iteration: ",x+1,"/",num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ",len(verts),len(faces))
        A=verts[faces[0][0]]
        B=verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A,B)*Rt/1000.)))
        print("   - Selecting faces on bounds ",bounds)

        lfaces=[]
        for face in faces:
            #print "#face: ",face
            lgeo=[]
            for j in [0,1,2]:
                (x,y,z)=verts[face[j]]
                (rlon,rlat,rh)=coordinates.xyz2geospheric(x,y,z)
                lgeo.append((rlon,rlat))
            
 #           triangle=Polygon.Polygon((lgeo))
 #           intersection = rectangle & triangle

            shp_triangle = shp_polygon(lgeo)
            shp_intersection = shp_rectangle.intersects(shp_triangle)

#            if not intersection:continue#print "No intersection"
            if not shp_intersection: continue  # print "No intersection"
            #else: print "Intersection"
            lfaces.append(face)
        
    faces=lfaces
    print("-- ",len(faces)," faces kept.")
    
    return (verts,faces)

    
    
    ftriangles=open('triangles.dat','w')
    for i in range(len(faces)):
        print("- face #",i,"/",len(faces))
        ftriangles.write(">\n")
        face=faces[i]
        for j in [0,1,2]:
            (x,y,z)=verts[face[j]]
            #print 'x,y,z ',x,y,z
            x=x*Rt
            y=y*Rt
            z=z*Rt
            (l,phi,he)=coordinates.xyz2geo(x,y,z)
            #print ("%10.5lf %10.5lf \n" % (math.degrees(l),math.degrees(phi)))
            ftriangles.write("%10.5lf %10.5lf \n" % (math.degrees(l),math.degrees(phi)))
        
    ftriangles.close()

    Rt=6371.0E3
    (verts,faces)=icosahedron()
    print("-- Number of vertices and faces of initial icosahedron ",len(verts),len(faces))
    print("-- Number of divisions ",num_subdivisions)
    print("-- Now doing subdivision...")
    for x in range(num_subdivisions):
        print("   - Division iteration: ",x+1,"/",num_subdivisions)
        verts, faces = subdivide(verts, faces)
        print("   - New Number of vertices and faces ",len(verts),len(faces))
        A=verts[faces[0][0]]
        B=verts[faces[0][1]]
        print(("   - New triangle edge distance : %8.3lf km" % (distance(A,B)*Rt/1000.)))
    return (verts,faces)
    
    
    
#    ftriangles=open('triangles.dat','w')
#    for i in range(len(faces)):
#        print "- face #",i,"/",len(faces)
#        ftriangles.write(">\n")
#        face=faces[i]
#        for j in [0,1,2]:
#            (x,y,z)=verts[face[j]]
#            #print 'x,y,z ',x,y,z
#            x=x*Rt
#            y=y*Rt
#            z=z*Rt
#            (l,phi,he)=Coordinates.xyz2geo(x,y,z)
#            #print ("%10.5lf %10.5lf \n" % (math.degrees(l),math.degrees(phi)))
#            ftriangles.write("%10.5lf %10.5lf \n" % (math.degrees(l),math.degrees(phi)))
#        
#    ftriangles.close()
    
