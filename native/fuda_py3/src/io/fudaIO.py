# fudaIO module: Input/output helper module for fuda.

# fuda [FUnction and Data Analysis tool] is a user interface to the
# fudalib module. It wraps all the functions in the fudalib module and
# supplement with some new functions.

import dataIO, fuda, sys
from types import *


######################################################################
#
# Exceptions.
#
######################################################################

class fudaIOError(Exception):
    def __init__(self, message):
        self.message = message


######################################################################
#
# Wrapper functions.
#
######################################################################

def open(*arguments,**keywords):
    return dataIO.open(*arguments,**keywords)

def close(*arguments,**keywords):
    return dataIO.close(*arguments,**keywords)

def info(*arguments,**keywords):
    return dataIO.info(*arguments,**keywords)

def set_pt_offset(*arguments,**keywords):
    return dataIO.set_pt_offset(*arguments,**keywords)

def set_dim_offset(*arguments,**keywords):
    return dataIO.set_dim_offset(*arguments,**keywords)

def get_pt_offset(*arguments,**keywords):
    return dataIO.get_pt_offset(*arguments,**keywords)

def get_dim_offset(*arguments,**keywords):
    return dataIO.get_dim_offset(*arguments,**keywords)

def get_dim(*arguments,**keywords):
    return dataIO.get_dim(*arguments,**keywords)

def get_size(*arguments,**keywords):
    return dataIO.get_size(*arguments,**keywords)

def read_pt(*arguments,**keywords):
    return dataIO.read_pt(*arguments,**keywords)

def read_mx(*arguments,**keywords):
    return dataIO.read_mx(*arguments,**keywords)

def map_default(*arguments,**keywords):
    return dataIO.map_default(*arguments,**keywords)

def map_linear(*arguments,**keywords):
    return dataIO.map_linear(*arguments,**keywords)

def map_point(*arguments,**keywords):
    return dataIO.map_point(*arguments,**keywords)

def map_p2x(*arguments,**keywords):
    return dataIO.map_p2x(*arguments,**keywords)

def map_x2p(*arguments,**keywords):
    return dataIO.map_x2p(*arguments,**keywords)


######################################################################
#
# Python extension functions functions.
#
######################################################################


def map_xs2ps(io_index,xseq):

    # Convert a user unit vector (sequence) to a point sequence and
    # return it.

    # Make point list.
    pseq=[]
    for idim in range(dataIO.get_dim(io_index)):
        pseq.append(map_x2p(io_index,idim,xseq[idim]))

    return pseq


def map_ps2xs(io_index,pseq):

    # Convert a point vector (sequence) to a user unit sequence and
    # return it.

    # Make user unit list.
    xseq=[]
    for idim in range(dataIO.get_dim(io_index)):
        xseq.append(map_p2x(io_index,idim,pseq[idim]))

    return xseq

def TupleType(key):
    return isinstance(key, (tuple, list))

def data_read_pt(io_index,coord,u):

    # Read a data point from the file with io_index as specified by
    # the coord sequence and define a fuda data point for the current
    # function with the coord translated to user units and with the
    # uncertainty u.

    # Get dimension of file:
    dim=dataIO.get_dim(io_index)

    # Check coord type and dimension.
    if TupleType(coord):# or type(coord)==ListType:
        pass
    else:
        raise fudaIOError('data_read_pt: invalid second argument')

    if not len(coord)==dim:
        raise fudaIOError('data_read_pt: invalid second argument')

    # Read the point in the file.
    data_value=dataIO.read_pt(io_index,coord)

    # Convert coord points to user units.
    ucoord=list(map_ps2xs(io_index,coord))

    # Append the data value and the uncertainty.
    ucoord.append(data_value)
    ucoord.append(u)

    # Add the data point to the current function in fuda.
    #apply(fuda.data,tuple(ucoord))
    fuda.data(*ucoord)


def data_read_center(io_index,center,radius,u,shape='r'):

    # This function reads a region of specified shape in an
    # N-dimensional data file centered at the point coordinates given
    # in the sequence center and with a radius in the N dimensions
    # given in the radius sequence. The radius sequence specifies a
    # number of points. E.g in a two-dimensional file, the call
    # data_read_center(io_index,(30,40),(5,6)) will read in
    # (2*5+1)*(2*6+1) data points and the centre point of the region
    # will be the point (30,40). The data points coordinates are
    # converted to user units and are added to the current fuda
    # function with an uncertainty of u.  The shape is a string of 'r'
    # and 'e' for rectangle and ellipsis, respectively. There should
    # be as many characters as there are dimensions in the file. If
    # there are less, the last character is repeated for the last
    # dimensions. If 'e' is specified it must be specified for at
    # least two dimensions to make sense.

    # Get dimension and size of file:
    dim=dataIO.get_dim(io_index)
    fsize=[]
    for idim in range(dim):
        fsize.append(dataIO.get_size(io_index,idim))

    # Get point offset.
    pt_offset=get_pt_offset(io_index)

    # Check center and radius type and dimension.
    if TupleType(center): # or type(center)==ListType:
        pass
    else:
        raise fudaIOError('data_read_center: invalid second argument')

    if not len(center)==dim:
        raise fudaIOError('data_read_center: invalid second argument')

    if TupleType(radius): # or type(radius)==ListType:
        pass
    else:
        raise fudaIOError('data_read_center: invalid 3rd argument')

    if not len(radius)==dim:
        raise fudaIOError('data_read_center: invalid second argument')

    # Check the shape.
    if not type(shape)==str:
        raise fudaIOError('data_read_center: shape must be a string')
    if len(shape)>dim:
        raise fudaIOError('data_read_center: shape string too long')

    # Eventyally, repeat last shape code.
    while len(shape)<dim:
        l=len(shape)
        shape=shape+shape[l-1]

    # Setup elipse index list containing the dimensions to define the
    # ellipses. Also, check shape string.
    ellipse=[]
    for i in range(dim):
        if shape[i]=='e':
            ellipse.append(i)
        elif shape[i]=='r':
            pass
        else:
            raise fudaIOError('data_read_center: invalid shape string')

    # Set ellipse dimension.
    ellipse_dim=len(ellipse)

    # We can not have a single dimension make up an ellipse.
    if ellipse_dim==1:
        raise fudaIOError('data_read_center: shape can not contain single e')

    # Convert center to the beginning of the rectangle.
    first=[]
    for idim in range(dim):
        first.append(center[idim]-radius[idim])

    # Calculate the size of the rectangle.
    size=[]
    for idim in range(dim):
        size.append(2*radius[idim]+1)

    # Initialize the relative and absolute coordinates
    # (the absolute coord is first[]+rcoord[]).
    rcoord=[]
    coord=[]
    for idim in range(dim):
        rcoord.append(0)
        coord.append(0)

    # Initialize the fuda data entry which needs space for the coordinate,
    # the point value and the uncertainty.
    data_entry=[]
    for idim in range(dim+2):
        data_entry.append(0.0)

    # Loop over all points of all dimensions.
    finido=0
    while not finido:

        # Calculate absolute coordinate.
        for idim in range(dim):
            coord[idim]=first[idim]+rcoord[idim]

        # Are we within bounds in the file.
        withinbounds=1
        for idim in range(dim):
            if coord[idim]<pt_offset or coord[idim]>fsize[idim]-1+pt_offset:
                withinbounds=0

        # Are we within the shape.
        if ellipse_dim>0:
            sum=0.0
            for i in range(ellipse_dim):
                idim=ellipse[i]
                sum=sum+pow(float(coord[idim]-center[idim])/\
                            (float(radius[idim])+0.5),2)
            if sum>1.0:
                withinbounds=0
        else:
            pass

        # Read the data if appropriate.
        if withinbounds:

            # Read from file and register data point with current fuda
            # function.
            data_read_pt(io_index,tuple(coord),u)
            
        # Step rcoord to next point.
        for idim in range(dim):
            if rcoord[idim]<size[idim]-1:
                # We increment this dimension and bread out.
                rcoord[idim]=rcoord[idim]+1
                break
            else:
                if idim==dim-1:

                    # If final dimension, we are done.
                    finido=1
                else:
                    # reset this rcoord.
                    rcoord[idim]=0



def data_read_xcenter(io_index,xcenter,radius,u,shape='r'):

    # This function is the same as data_read_center, except the
    # function takes the center in user units.

    # Call data_read_center with xcenter converted to points.
    data_read_center(io_index,map_xs2ps(io_index,xcenter),radius,u,shape)


