"""Helpers for generating and executing nmrPipe processing scripts.

The processing GUI uses the ProcessFrame builders to render either the plain
processing script (nmrproc.test.com) or the LP variant (nmrprocLP.com).
This module keeps the high-level create/save/execute orchestration in one place
so the GUI can preview the full script text without relying on a partially
written file on disk.
"""

from __future__ import annotations

import logging
import subprocess
import os
import shlex
from pathlib import Path
from tempfile import TemporaryDirectory
from spinDecon.processing.script_context import ProcessingScriptState, PROCESSING_CONTROL_NAMES
import nmrglue as ng
from typing import Optional, Tuple
import copy
import numpy



def _report_processing_error(message):
    """Report a processing validation error without depending on GUI infrastructure."""
    logging.error('%s', message)


def _path_escape(indir):
    fields = indir.split(' ')
    logging.info(fields)
    final_string = ''
    for x in range(len(fields) - 1):
        if fields[x][-1] != '\\':
            final_string += fields[x] + '\\ '
        else:
            final_string += fields[x]

    final_string += fields[-1]
    return final_string


def MakeProj4D(infile):
    #make projections
    indir,filename = os.path.split(infile)
    indir = _path_escape(indir)
    infile=  _path_escape(infile)
    with tempfile.TemporaryDirectory(prefix='decon_proj_') as tmpdir:
        script_path = os.path.join(tmpdir, 'proj.com')
        outy=open(script_path,'w')
        outy.write('echo \'Making projections\'\\n')
        outy.write('if( -e ft) rm -rf ft\\n')
        outy.write('nmrPipe -in %s  \\\\n' % infile)
        outy.write('|pipe2xyz -out ft/test%03d%03d.ft4 -x\\n')
        outy.write('rm *.dat\\n')
        outy.write('proj4D.tcl -in ft/test%03d%03d.ft4\\n')
        outy.write('rm -rf %s/projections\\n' % (indir))
        outy.write('mkdir %s/projections \\n' % (indir))
        outy.write('mv *.dat %s/projections \\n' % (indir))
        #outy.write('mv %s %s\\n' % (inf,indir+'/test.ft4'))
        outy.write('if( -e ft) rm -rf ft\\n')
        outy.write('echo Done.\\n')
        outy.close()
        os.system('csh %s' % _path_escape(script_path))
    return 1



def _projection_mode(mode):
    """Return the supported projection reduction mode, defaulting to sum."""
    value = str(mode or 'sum').strip().lower()
    return value if value in ('sum', 'skyline') else 'sum'


def _project_data(data, axis, mode='sum'):
    """Reduce *data* by summation or signed maximum-absolute skyline.

    Skyline preserves the sign of the point having the largest absolute
    intensity across all projected-out axes.
    """
    mode = _projection_mode(mode)
    if mode == 'sum':
        return numpy.sum(data, axis=axis)

    axes = (axis,) if isinstance(axis, int) else tuple(axis)
    axes = tuple(a if a >= 0 else data.ndim + a for a in axes)
    keep_axes = tuple(a for a in range(data.ndim) if a not in axes)
    ordered = keep_axes + axes
    work = numpy.transpose(data, ordered)
    keep_shape = tuple(data.shape[a] for a in keep_axes)
    reduce_size = int(numpy.prod([data.shape[a] for a in axes]))
    work = work.reshape(keep_shape + (reduce_size,))
    indices = numpy.argmax(numpy.abs(work), axis=-1)
    return numpy.take_along_axis(work, indices[..., numpy.newaxis], axis=-1)[..., 0]

def nmrglue_project3D(dic, data, dim=0, folder='projections', suffix='', projection_type='sum', progress=None):
    
    if int(dic['FDDIMCOUNT']) != 3:
        _report_processing_error('Error: 3D data expected')
        return
    # exit()
    dimorders = dic['FDDIMORDER']
    # remove 4 from dimorders
    dimorders = [x for x in dimorders if x != 4.0]
    H_size = dic['FDSIZE']
    i = 2.0
    for dim in dimorders:
        dim = dim
        newdic = dic.copy()
        newdic['FDDIMCOUNT'] = '2' # set dimension count to 2

        # move the current dimension to the third position
        dimorders2 = copy.deepcopy(dimorders)
        dimorders2.remove(dim)
        dimorders2.append(dim)
        dimorders2.append(4.0)
        # print(dimorders2)
        newdic['FDDIMORDER'] = dimorders2 
        newdic['FDDIMORDER1'] = dimorders2[0]
        newdic['FDDIMORDER2'] = dimorders2[1]
        newdic['FDDIMORDER3'] = dimorders2[2]
        newdic['FDDIMORDER4'] = dimorders2[3]

        
        # for key in newdic.keys():
        #     if 'FDF'+str(int(dim)) in key:
        #         print(key, dic[key])
            ## we want to change the order of the dimensions by changing the FDDIMORDER of the newdic.  As we are projecting the current directory, we want to move that dimension to position three, and the other two to 1 and 2.  we have to leave 4 at the end of the list (position 4).
        dim1 = dimorders2[0]
        label1 = newdic['FDF'+str(int(dim1))+'LABEL']
        dim2 = dimorders2[1]
        label2 = newdic['FDF'+str(int(dim2))+'LABEL']
        newdic['FDSPECNUM'] = newdic['FDF'+str(int(dim2))+'FTSIZE']
        newdic['FDSIZE'] = newdic['FDF'+str(int(dim1))+'FTSIZE']

        newdic['FDF'+str(int(dim))+'FTSIZE'] = 0
        
        if dim2 == 2.0:
            newdic['FDSPECNUM'] = H_size
        if dim1 == 2.0:
            newdic['FDSIZE'] = H_size
        
        # actual_dim_location = dimorders.index(i)
        actual_dim_location = int(i)
        # print(data.shape, actual_dim_location)
        ## find either max or min of the data along the current dimension

        data_proj1 = numpy.max(data, axis=actual_dim_location)
        data_proj2 = numpy.max(-data, axis=actual_dim_location)

        # Sum the data along the projected-out axis so the saved projection
        # preserves integrated intensity rather than peak height.
        data_proj = _project_data(data, axis=actual_dim_location, mode=projection_type)
        

        # print(data_proj.shape)
        
        # print(label1, label2)
        # print(dim1, dim2)
        # print(dic['FDF1FTSIZE'], dic['FDF2FTSIZE'], dic['FDF3FTSIZE'])

        # print(newdic['FDSIZE'], newdic['FDSPECNUM'])

        outfile = os.path.join(folder, '%s.%s%s.dat' % (label1, label2, suffix))
        ng.pipe.write(outfile, newdic, data_proj, overwrite=True)
        if progress:
            progress('  Created 2D projection: %s x %s\n' % (label1, label2))
        # print()
        i-=1.0

def nmrglue_project3p_spectral(dic, data, folder='projections', suffix='', progress=None):
    """Create the single direct/indirect spectral projection for ``3p`` data.

    The pseudo/real axis is summed (never skyline projected).  The saved array
    is normalised to the processProjectionFrame convention: rows are the
    direct spectral dimension (F2; Y axis) and columns are the indirect
    spectral dimension (F1; X axis).  The filename follows the same convention,
    e.g. ``H1.C13.dat``.
    """
    if int(float(dic.get('FDDIMCOUNT', 0))) != 3:
        _report_processing_error('Error: pseudo-3D projection requires 3D data')
        return None

    arr = numpy.asarray(data)
    if arr.ndim != 3:
        _report_processing_error('Error: pseudo-3D projection requires a 3D data array')
        return None

    dimorders = [float(x) for x in dic.get('FDDIMORDER', []) if float(x) != 4.0]
    if len(dimorders) != 3:
        _report_processing_error('Error: could not determine pseudo-3D NMRPipe dimension order')
        return None

    pseudo_dims = []
    for d in dimorders:
        try:
            if float(dic.get('FDF%dFTFLAG' % int(d), 0.0)) == 0.0:
                pseudo_dims.append(d)
        except (TypeError, ValueError):
            pass

    if len(pseudo_dims) != 1:
        pseudo_labels = ('time_', 'ncyc', 'gzlvl', 'id', 'relax', 'delay')
        labelled = []
        for d in dimorders:
            label = str(dic.get('FDF%dLABEL' % int(d), '')).strip().lower()
            if any(token in label for token in pseudo_labels):
                labelled.append(d)
        if len(labelled) == 1:
            pseudo_dims = labelled

    if len(pseudo_dims) != 1:
        labels = [str(dic.get('FDF%dLABEL' % int(d), d)) for d in dimorders]
        _report_processing_error('Error: could not uniquely identify real/pseudo axis for projection: %s' % labels)
        return None

    pseudo_dim = pseudo_dims[0]
    pseudo_position = dimorders.index(pseudo_dim)
    pseudo_axis = arr.ndim - 1 - pseudo_position

    data_proj = numpy.sum(arr, axis=pseudo_axis)
    if data_proj.ndim != 2:
        _report_processing_error('Error: pseudo-3D spectral projection is not two-dimensional')
        return None

    # In this acquisition type the two retained spectral dimensions are the
    # conventional NMRPipe F2 (direct) and F1 (indirect) dimensions.  Work out
    # where each one sits in the projected NumPy array from FDDIMORDER rather
    # than assuming the pseudo dimension occupied a particular NumPy axis.
    direct_dim = 2.0 if 2.0 in dimorders and 2.0 != pseudo_dim else None
    indirect_dim = 1.0 if 1.0 in dimorders and 1.0 != pseudo_dim else None
    spectral_dims = [d for d in dimorders if d != pseudo_dim]
    if direct_dim is None or indirect_dim is None:
        # Conservative fallback for unusual headers: FDDIMORDER[0] is the
        # direct/last NumPy dimension and the other retained dimension is indirect.
        direct_dim = spectral_dims[0]
        indirect_dim = spectral_dims[1]

    remaining_positions = [i for i, d in enumerate(dimorders) if d != pseudo_dim]
    remaining_numpy_axes = [arr.ndim - 1 - i for i in remaining_positions]
    # numpy.sum removes pseudo_axis; map original surviving axes into the 2D result.
    original_to_projected = {}
    out_axis = 0
    for original_axis in range(arr.ndim):
        if original_axis == pseudo_axis:
            continue
        original_to_projected[original_axis] = out_axis
        out_axis += 1
    dim_to_projected_axis = {}
    for d, original_axis in zip(spectral_dims, remaining_numpy_axes):
        dim_to_projected_axis[d] = original_to_projected[original_axis]

    direct_axis = dim_to_projected_axis[direct_dim]
    indirect_axis = dim_to_projected_axis[indirect_dim]
    # Match the conventional NMRPipe 2D projection layout used by the normal
    # 3D projector: rows/Y are the indirect spectral dimension (F1), while
    # columns/X are the direct spectral dimension (F2).
    if (indirect_axis, direct_axis) != (0, 1):
        data_proj = numpy.moveaxis(data_proj, (indirect_axis, direct_axis), (0, 1))

    # Build a self-consistent conventional 2D NMRPipe header.  nmrglue maps
    # NumPy axis 0 -> F1 (rows/Y) and axis 1 -> F2 (columns/X), so F1 must
    # carry the indirect metadata and F2 the direct metadata.
    source_dic = dic.copy()
    newdic = dic.copy()

    direct_label = str(source_dic.get('FDF%dLABEL' % int(direct_dim), 'Direct')).strip()
    indirect_label = str(source_dic.get('FDF%dLABEL' % int(indirect_dim), 'Indirect')).strip()

    def _copy_fdf_metadata(src_dim, dst_dim):
        src_prefix = 'FDF%d' % int(src_dim)
        dst_prefix = 'FDF%d' % int(dst_dim)
        snapshot = {k: v for k, v in source_dic.items() if k.startswith(src_prefix)}
        # Keep the destination dimension's complete NMRPipe field set intact.
        # nmrglue.dic2fdata() expects a number of mandatory FDF1/FDF2 keys
        # (for example APODDF) even when the corresponding source dimension
        # did not explicitly carry that key.  Earlier code deleted the whole
        # destination prefix before copying and could therefore create an
        # incomplete dictionary.  Overlay the source values instead; this
        # preserves required defaults while still remapping all metadata that
        # is actually present on the source spectral dimension.
        for key, value in snapshot.items():
            newdic[dst_prefix + key[len(src_prefix):]] = value

    # Physical output axis 0 (rows/Y) is indirect -> NMRPipe F1.
    # Physical output axis 1 (columns/X) is direct -> NMRPipe F2.
    _copy_fdf_metadata(indirect_dim, 1.0)
    _copy_fdf_metadata(direct_dim, 2.0)

    new_order = [2.0, 1.0, 3.0, 4.0]
    newdic['FDDIMCOUNT'] = 2.0
    # Conventional 2D projection files are not transposed.
    newdic['FDTRANSPOSED'] = 0.0
    newdic['FDDIMORDER'] = new_order
    for idx, value in enumerate(new_order, start=1):
        newdic['FDDIMORDER%d' % idx] = value

    # For a 2D NMRPipe file FDSIZE is the F2/column count and FDSPECNUM is
    # the F1/row count.  Keep all size metadata consistent with data_proj.
    newdic['FDSIZE'] = int(data_proj.shape[1])
    newdic['FDSPECNUM'] = int(data_proj.shape[0])
    newdic['FDF1FTSIZE'] = int(data_proj.shape[0])
    newdic['FDF2FTSIZE'] = int(data_proj.shape[1])
    newdic['FDF1TDSIZE'] = int(data_proj.shape[0])
    newdic['FDF2TDSIZE'] = int(data_proj.shape[1])
    # The third dimension is no longer part of this 2D file.
    newdic['FDF3FTSIZE'] = 0.0
    outfile = os.path.join(folder, '%s.%s%s.dat' % (direct_label, indirect_label, suffix))
    os.makedirs(folder, exist_ok=True)
    ng.pipe.write(outfile, newdic, data_proj, overwrite=True)

    check_dic, check_data = ng.pipe.read(outfile)
    if tuple(numpy.asarray(check_data).shape) != tuple(data_proj.shape):
        raise RuntimeError(
            'Pseudo-3D projection header/data mismatch for %s: wrote %s, read %s'
            % (outfile, data_proj.shape, numpy.asarray(check_data).shape)
        )

    logging.info(
        'Created pseudo-3D spectral sum projection %s by summing axis %d (%s); '
        'orientation rows=%s columns=%s shape=%s',
        outfile, pseudo_axis, str(newdic.get('FDF%dLABEL' % int(pseudo_dim), pseudo_dim)),
        indirect_label, direct_label, data_proj.shape,
    )
    if progress:
        progress('  Created pseudo-3D spectral sum projection: %s (Y) x %s (X)\n' %
                 (indirect_label, direct_label))
    return outfile


def MakeProj3P(infile, folder='projections', clean=True, suffix='', progress=None):
    """Create only the two-spectral-dimension projection for ``3p`` data."""
    logging.info('Generating pseudo-3D spectral projection from %s', infile)
    dic, data = ng.pipe.read(infile)
    if not folder:
        indir = os.path.dirname(infile)
        folder = os.path.join(indir, 'projections') if indir else 'projections'
    os.makedirs(folder, exist_ok=True)
    if clean:
        for filename in os.listdir(folder):
            path = os.path.join(folder, filename)
            if os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass
    # First create the single canonical 2D spectral projection.  For 3p data
    # this is the only meaningful 2D projection: the pseudo/real axis has
    # already been summed out by nmrglue_project3p_spectral().
    projection_file = nmrglue_project3p_spectral(
        dic, data, folder=folder, suffix=suffix, progress=progress
    )
    if not projection_file:
        return projection_file

    # The rest of the GUI expects the direct-dimension 1D projection (for
    # example H1.dat) in spec/projections1D.  Generate both spectral 1D
    # projections from the *canonical 2D projection*, rather than directly
    # from the pseudo-3D cube.  This keeps their axes/header conventions
    # identical to H1.C13.dat and deliberately excludes the pseudo axis.
    folder1D = os.path.join(os.path.dirname(folder), 'projections1D')
    os.makedirs(folder1D, exist_ok=True)
    if clean:
        for filename in os.listdir(folder1D):
            path = os.path.join(folder1D, filename)
            if os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

    proj_dic, proj_data = ng.pipe.read(projection_file)
    nmrglue_project2D_1D(
        proj_dic, proj_data, folder=folder1D, projection_type='sum', progress=progress
    )
    return projection_file


def nmrglue_project3D_1D(dic, data, dim=0, folder='projections', projection_type='sum', progress=None):
    
    if int(dic['FDDIMCOUNT']) != 3:
        _report_processing_error('Error: 3D data expected')
        return
    # exit()
    dimorders = dic['FDDIMORDER']
    # remove 4 from dimorders
    dimorders = [x for x in dimorders if x != 4.0]
    H_size = dic['FDSIZE']
    i = 2.0
    for dim in dimorders:
        dim = dim
        newdic = dic.copy()
        newdic['FDDIMCOUNT'] = '1' # set dimension count to 1

        dimorders2 = copy.deepcopy(dimorders)
        dimorders2.remove(dim)
        dimorders2.append(dim)
        dimorders2.append(4.0)
  
        dim1 = dimorders2[0]
        label1 = newdic['FDF'+str(int(dim1))+'LABEL']
        dim2 = dimorders2[1]
        label2 = newdic['FDF'+str(int(dim2))+'LABEL']
 
        dim3 = dim
        label3 = newdic['FDF'+str(int(dim3))+'LABEL']
        newdic['FDSIZE'] = newdic['FDF'+str(int(dim3))+'FTSIZE']
        if dim3 == 2.0:
            newdic['FDSIZE'] = H_size
        newdic['FDSPECNUM'] = 0

        dimorders2=dim3,dim1,dim2,4.0

        newdic['FDDIMORDER'] = dimorders2 
        newdic['FDDIMORDER1'] = dimorders2[0]
        newdic['FDDIMORDER2'] = dimorders2[1]
        newdic['FDDIMORDER3'] = dimorders2[2]
        newdic['FDDIMORDER4'] = dimorders2[3]

        #take the two that we don't want...
        actual_dim_location = (int(i-1)%3,int(i+1)%3)
 
        data_proj = _project_data(data, axis=actual_dim_location, mode=projection_type)
        ng.pipe.write(folder+'/%s.dat' % (label3), newdic, data_proj, overwrite=True)
        if progress:
            progress('  Created 1D projection: %s\n' % label3)
        # print()
        i-=1.0

def nmrglue_project2D_1D(dic, data, dim=0, folder='projections1D', projection_type='sum', progress=None):
    
    if int(dic['FDDIMCOUNT']) != 2:
        _report_processing_error('Error: 2D data expected')
        return
    # exit()
    dimorders = dic['FDDIMORDER']
    # remove 4 from dimorders
    dimorders = [x for x in dimorders if x != 4.0]

    H_size = dic['FDSIZE']
    i = 1.0
    for dim in dimorders:
        dim = dim
        newdic = dic.copy()
        newdic['FDDIMCOUNT'] = '1' # set dimension count to 1

        dimorders2 = copy.deepcopy(dimorders)
        dimorders2.remove(dim)
        dimorders2.append(dim)
        dimorders2.append(4.0)
  
        dim1 = dimorders2[0]
        label1 = newdic['FDF'+str(int(dim1))+'LABEL']
        dim2 = dimorders2[1]
        label2 = newdic['FDF'+str(int(dim2))+'LABEL']
 
        dim3 = dim
        label3 = newdic['FDF'+str(int(dim3))+'LABEL']
        newdic['FDSIZE'] = newdic['FDF'+str(int(dim3))+'FTSIZE']     
        if dim3 == 2.0:
            newdic['FDSIZE'] = dic['FDSPECNUM']
        newdic['FDSPECNUM'] = 0

        dimorders2=dim3,dim1,dim2,4.0
        newdic['FDDIMORDER'] = dimorders2 
        newdic['FDDIMORDER1'] = dimorders2[0]
        newdic['FDDIMORDER2'] = dimorders2[1]
        newdic['FDDIMORDER3'] = dimorders2[2]
        newdic['FDDIMORDER4'] = dimorders2[3]

        #take the two that we don't want...
        actual_dim_location = (int(i+1)%2,)
 
        data_proj = _project_data(data, axis=actual_dim_location, mode=projection_type)

        # A 1D NMRPipe file must advertise the number of points actually
        # written.  In 2D spectra the direct dimension is the final numpy
        # axis (data.shape[-1]).  Older code special-cased FDF2 using
        # FDSPECNUM, which is the indirect dimension size for normal 2D
        # data; this produced e.g. a 1023-point 1H projection with a 2048
        # point header.  Use the projected trace itself as the authority.
        newdic['FDSIZE'] = int(numpy.asarray(data_proj).size)
        if dim3 == 2.0:
            direct_size = int(data.shape[-1])
            if newdic['FDSIZE'] != direct_size:
                logging.warning(
                    'Direct-dimension projection size mismatch: projection=%s, direct=%s',
                    newdic['FDSIZE'], direct_size,
                )

        ng.pipe.write(folder+'/%s.dat' % (label3), newdic, data_proj, overwrite=True)
        if progress:
            progress('  Created 1D projection: %s\n' % label3)
        # print()
        i-=1.0
        if(i<0):
            break

def MakeProj3D(infile, folder='projections', OneD=False, clean=True, suffix='', projection_type='sum', progress=None):
    #make projections
    logging.info('Generating 3D projections from %s', infile)
    indir,filename = os.path.split(infile)
    indir = _path_escape(indir)
    infile=  _path_escape(infile)
    dic, data = ng.pipe.read(infile)
    ## project 3D data onto 3 2D planes using nmrglue rather than nmrpipe, which is now redundant
    # 1. project onto XY plane
    if not folder:
        folder = os.path.join(indir, 'projections') if indir else 'projections'
    os.makedirs(folder, exist_ok=True)
    if clean:
        for filename in os.listdir(folder):
            path = os.path.join(folder, filename)
            if os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass
    proj = nmrglue_project3D(dic, data, dim=0, folder=folder, suffix=suffix, projection_type=projection_type, progress=progress)

    if(OneD==False):
        return
    folder1D=folder.replace('projections','projections1D')
    os.makedirs(folder1D, exist_ok=True)
    for filename in os.listdir(folder1D):
        path = os.path.join(folder1D, filename)
        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass
    proj = nmrglue_project3D_1D(dic, data, dim=0, folder=folder1D, projection_type=projection_type, progress=progress)

    # outy=open('raw/proj.com','w')
    # outy.write('echo \'Making projections\'\n')
    # outy.write('if( -e ft) rm -rf ft\n')
    # outy.write('nmrPipe -in %s  \\\n' % infile)
    # outy.write('|pipe2xyz -out ft/test%03d.ft3 -x\n')
    # outy.write('rm *.dat\n')
    # outy.write('proj3D.tcl -in ft/test%03d.ft3\n')
    # outy.write('rm -rf %s/projections\n' % (indir))
    # outy.write('mkdir %s/projections \n' % (indir))
    # outy.write('mv *.dat %s/projections \n' % (indir))
    # #outy.write('mv %s %s\n' % (inf,self._spec_output_dir()+'/test.ft3'))
    # outy.write('if( -e ft) rm -rf ft\n')
    # outy.write('echo Done.\n')
    # outy.close()
    # os.system('csh raw/proj.com')

    #outy.write('rm -rf {}/projections1D\n'.format(self._spec_output_dir()))
    #outy.write('mkdir {}/projections1D\n'.format(self._spec_output_dir()))
    #    outy.write('cd {}/projections\n'.format(self._spec_output_dir()))
    #    outy.write('foreach file (*.dat)\n')
    #    outy.write('    cp $file ../projections1D\n')
    #    outy.write('end\n')
    #    outy.write('cd ../projections1D\n')
    #    outy.write('foreach f (*.dat)\n')
    #    outy.write('mv $f $f:r.ft2\n')
    #    outy.write('end\n')
    #    outy.write('foreach file (*.ft2)\n')
    #    outy.write('    proj2D.tcl -in $file\n')
    #    outy.write('end\n')
    #    outy.write('rm *.ft2\n')
    #    outy.write('cd ../..\n')

    return 1

def MakeProj2D(infile, projection_type='sum', progress=None):
    #make projections
    logging.info('Generating 2D projections from %s', infile)
    indir,filename = os.path.split(infile)
    indir = _path_escape(indir)
    infile=  _path_escape(infile)
    folder=os.path.join(indir,'projections')
    dic, data = ng.pipe.read(infile)
    ## project 3D data onto 3 2D planes using nmrglue rather than nmrpipe, which is now redundant
    # 1. project onto XY plane
    os.makedirs(folder, exist_ok=True)
    for filename in os.listdir(folder):
        path = os.path.join(folder, filename)
        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass
    proj = nmrglue_project2D_1D(dic, data, dim=0, folder=folder, projection_type=projection_type, progress=progress)
    return 1

def script_path_for(frame, lp: str = 'n') -> str:
    filename = 'nmrprocMDDNMR.com' if lp == 'm' else ('nmrprocLP.com' if lp == 'y' else 'nmrproc.test.com')
    try:
        base_dir = frame._spec_output_dir()
    except Exception:
        base_dir = ''
    if not base_dir:
        base_dir = './spec'
    return os.path.join(base_dir, filename)


def pipefile_for(frame) -> str:
    if frame.spectral_dim_count == 1:
        return 'test.ft'
    if frame.spectral_dim_count == 2:
        return 'test.ft2'
    if frame.spectral_dim_count == 3:
        return 'test.ft3'
    return 'test.ft4'




def _widget_value(widget, default=None):
    """Accept either a legacy wx control or a plain script-state value."""
    if widget is None:
        return default
    method = getattr(widget, 'GetValue', None)
    if callable(method):
        try:
            return method()
        except Exception:
            return default
    return widget


def _checked(widget) -> bool:
    """Accept either a legacy checkbox or a plain boolean value."""
    method = getattr(widget, 'IsChecked', None)
    if callable(method):
        try:
            return bool(method())
        except Exception:
            return False
    return bool(widget)


def get_xmin_xmax(frame):
    xmin_s = []
    xmax_s = []
    try:
        xmin_s = frame.xminBox.GetValue().split()
    except Exception:
        pass
    try:
        xmax_s = frame.xmaxBox.GetValue().split()
    except Exception:
        pass
    if len(xmin_s) > 0:
        if xmin_s[0] == '*':
            xmin = '*'
        else:
            xmin = float(xmin_s[0])
    else:
        xmin = '*'
    if len(xmax_s) > 0:
        if xmax_s[0] == '*':
            xmax = '*'
        else:
            xmax = float(xmax_s[0])
    else:
        xmax = '*'
    return xmin, xmax


def do_extract_f1(frame, outy):
    xmin, xmax = get_xmin_xmax(frame)
    if xmin == '*' and xmax == '*':
        return
    outy.write('| nmrPipe  -fn EXT ')
    if xmin != '*':
        outy.write('-x1 %sppm ' % xmin)
    else:
        outy.write('-x1 100%s ' % '%')
    if xmax != '*':
        outy.write('-xn %sppm ' % xmax)
    else:
        outy.write('-xn 0%s ' % '%')
    outy.write(' -sw \\\n')


def do_smile(frame, T, ax, p0, p1, f180, lp, wBox, w2Val, wVal, outy):
    if _checked(lp):
        T = 0
    p0v = float(_widget_value(p0, 0.0))
    p1v = float(_widget_value(p1, 0.0))
    if _checked(f180):
        p0v += -90
        p1v += 180

    win = _widget_value(wBox, 'SP')
    if win == 'SP':
        Q1 = float(_widget_value(wVal, 0.0))
        Q2 = 0.99
        Q3 = float(_widget_value(w2Val, 0.0))
    else:
        Q1 = float(_widget_value(wVal, 0.0))
        Q2 = float(_widget_value(w2Val, 0.0))
        Q3 = 0.0

    outy.write(' -%sT %f \\\n' % (ax, T))
    outy.write(' -%sP0 %f    -%sP1 %f \\\n' % (ax, p0v, ax, p1v))
    outy.write(' -%szf 2 \\\n' % ax)


def do_flags(frame, wbox, w2Val, w1Val, firstPoint, f180, lp, lp_glob, p0, p1, ft, outy):
    win = _widget_value(wbox, 'GM')
    win1Val = float(_widget_value(w1Val, 0.0))
    win2Val = float(_widget_value(w2Val, 0.0))
    c = float(_widget_value(firstPoint, 0.5))

    try:
        if _checked(f180):
            c = 1.0
    except Exception:
        pass

    try:
        if _checked(lp) and lp_glob == 'y':
            c = 1.0
    except Exception:
        pass

    p0v = float(_widget_value(p0, 0.0))
    p1v = float(_widget_value(p1, 0.0))

    try:
        if _checked(f180):
            p0v += -90
            p1v += 180
    except Exception:
        pass

    ft_choices = getattr(frame, 'ftlisty', ['Auto', 'Neg', 'Alt', 'AltNeg', 'Real'])
    try:
        selection_method = getattr(ft, 'GetSelection', None)
        selection = selection_method() if callable(selection_method) else int(ft or 0)
        ftsel = ft_choices[selection]
    except Exception:
        ftsel = 'Auto'
    if ftsel == 'Auto':
        ftf = '-auto'
    elif ftsel == 'Neg':
        ftf = '-neg'
    elif ftsel == 'Alt':
        ftf = '-alt'
    elif ftsel == 'AltNeg':
        ftf = '-alt -neg'
    elif ftsel == 'Real':
        ftf = '-real'
    else:
        ftf = ' '

    if win == 'SP':
        outy.write('| nmrPipe  -fn %s -off %f -end 0.99 -pow %f -c %f    \\\n' % (win, win1Val, win2Val, c))
    elif win == 'EM':
        outy.write('| nmrPipe  -fn %s -lb %f -c %f    \\\n' % (win, win1Val, c))
    else:
        outy.write('| nmrPipe  -fn %s -g1 %f -g2 %f -g3 0.0 -c %f    \\\n' % (win, win1Val, win2Val, c))
    if frame.spectral_dim_count == 4:
        outy.write('| nmrPipe  -fn ZF -auto                               \\\n')
    else:
        outy.write('| nmrPipe  -fn ZF -zf 2                               \\\n')
    outy.write('| nmrPipe  -fn FT %s                               \\\n' % ftf)
    outy.write('| nmrPipe  -fn PS -p0 %s -p1 %s -di              \\\n' % (p0v, p1v))


def make_1d_proj_3d(frame, outy):
    outy.write('rm -rf {}/projections1D\n'.format(frame._spec_output_dir()))
    outy.write('mkdir {}/projections1D\n'.format(frame._spec_output_dir()))
    outy.write('cd {}/projections\n'.format(frame._spec_output_dir()))
    outy.write('foreach file (*.dat)\n')
    outy.write('    cp $file ../projections1D\n')
    outy.write('end\n')
    outy.write('cd ../projections1D\n')
    outy.write('foreach f (*.dat)\n')
    outy.write('mv $f $f:r.ft2\n')
    outy.write('end\n')
    outy.write('foreach file (*.ft2)\n')
    outy.write('    proj2D.tcl -in $file\n')
    outy.write('end\n')
    outy.write('rm *.ft2\n')
    outy.write('cd ../..\n')


def make_proj_3dp(frame, inf, outy):
    outy.write('echo \'Making projections\'\n')
    outy.write('rm *.dat\n')
    outy.write('proj3D.tcl -in %s/test.ft2\n' % frame._spec_output_dir())
    outy.write('rm -rf %s/projections\n' % frame._spec_output_dir())
    outy.write('mkdir %s/projections \n' % frame._spec_output_dir())
    outy.write('mv *.dat %s/projections \n' % frame._spec_output_dir())
    outy.write('if( -e ft) rm -rf ft\n')
    outy.write('echo Done.\n')


def make_projy_3d(frame, inf, outy):
    outy.write('echo \'Making projections\'\n')
    outy.write('if( -e ft) rm -rf ft\n')
    outy.write('nmrPipe -in %s  \\\n' % inf)
    outy.write('|pipe2xyz -out ft/test%03d.ft3 -x\n')
    outy.write('rm *.dat\n')
    outy.write('proj3D.tcl -in ft/test%03d.ft3\n')
    outy.write('rm -rf %s/projections\n' % frame._spec_output_dir())
    outy.write('mkdir %s/projections \n' % frame._spec_output_dir())
    outy.write('mv *.dat %s/projections \n' % frame._spec_output_dir())


def make_proj_4d(frame, inf, outy):
    outy.write('echo \'Making projections\'\n')
    outy.write('if( -e ft) rm -rf ft\n')
    outy.write('nmrPipe -in %s  \\\n' % inf)
    outy.write('|pipe2xyz -out ft/test%03d%03d.ft4 -x\n')
    outy.write('rm *.dat\n')
    outy.write('proj4D.tcl -in ft/test%03d%03d.ft4\n')
    outy.write('rm -rf %s/projections\n' % frame._spec_output_dir())
    outy.write('mkdir %s/projections \n' % frame._spec_output_dir())
    outy.write('mv *.dat %s/projections \n' % frame._spec_output_dir())
    outy.write('mv %s %s\n' % (inf, frame._spec_output_dir() + '/test.ft4'))
    outy.write('if( -e ft) rm -rf ft\n')
    outy.write('echo Done.\n')
def make_proc_script_1d_state(process, state: ProcessingScriptState, outfile, lp='n'):
    logging.info('MakeProcScript1D')

    outy=open(outfile,'w')
    outy.write('#!/bin/csh -f\n')
    outy.write('%s\n' % ('set ft4trec='+process._spec_output_dir()+'/test.fid'))

    outy.write('echo \'Processing 1D\'\n')
    outy.write('showhdr $ft4trec\n')

    outy.write('echo Processing XY Dimensions Without LP ......\n')
    outy.write('if( -e %s/test.ft2) rm -rf %s/test.ft2\n' % (process._spec_output_dir(),process._spec_output_dir()))
    outy.write('xyz2pipe -in $ft4trec -verb -x \\\n')
    outy.write('# Process X \\\n')
    if(state.checked('cb_baseSol')):
        outy.write('| nmrPipe -fn SOL \\\n')

    do_flags(process,state.control('windowBox0'),state.control('win2Val0'),state.control('win3Val0'),state.control('firstPoint0'),False,False,'n',state.control('p0'),state.control('p1'),state.control('cb_ft0'),outy)

    do_extract_f1(process, outy)

    if(state.checked('cb_baseLin')):
        outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
    if(state.checked('cb_basepol')):
        outy.write('| nmrPipe -fn POLY -auto\\\n')

    outy.write('| pipe2xyz -out %s -x -verb -ov\n' % (process._spec_output_dir()+'/test.ft',))
    outy.close()

    return


def make_proc_script_1d(frame, outfile, lp='n'):
    """Legacy frame-shaped entry point retained during staged migration."""
    state = ProcessingScriptState.capture(frame, {}, PROCESSING_CONTROL_NAMES)
    return make_proc_script_1d_state(frame, state, outfile, lp=lp)


def make_proc_script_1d_slice(frame, outfile, lp='n'):
    logging.info('MakeProcScript1DSlice')

    outy=open(outfile,'w')
    outy.write('#!/bin/csh -f\n')
    outy.write('%s\n' % ('set ft4trec='+frame._spec_output_dir()+'/slice.fid'))

    outy.write('echo \'Processing 1D slice\'\n')
    outy.write('showhdr $ft4trec\n')

    outy.write('echo Processing XY Dimensions Without LP ......\n')
    outy.write('if( -e %s/slice.ft1) rm -rf %s/slice.ft1\n' % (frame._spec_output_dir(),frame._spec_output_dir()))
    outy.write('xyz2pipe -in $ft4trec -verb -x \\\n')

    # Apply the first-dimension apodisation settings exactly as configured in
    # the processing window. Prefer the live widgets, but fall back to the save
    # file/defaults if the processing frame is not currently open.
    if hasattr(frame, '_direct_apodization_settings'):
        win, win1Val, win2Val, c = frame._direct_apodization_settings()
    else:
        win = frame.windowBox0.GetValue() if hasattr(frame, 'windowBox0') else 'GM'
        try:
            win1Val = float(frame.win2Val0.GetValue()) if hasattr(frame, 'win2Val0') else 2.0
        except Exception:
            win1Val = 2.0
        try:
            win2Val = float(frame.win3Val0.GetValue()) if hasattr(frame, 'win3Val0') else 2.0
        except Exception:
            win2Val = 2.0
        try:
            c = float(frame.firstPoint0.GetValue()) if hasattr(frame, 'firstPoint0') else 0.5
        except Exception:
            c = 0.5
    try:
        if hasattr(frame, 'cb_baseSol') and frame.cb_baseSol.IsChecked():
            outy.write('| nmrPipe -fn SOL \\\n')
    except Exception:
        pass

    if win == 'SP':
        outy.write('| nmrPipe  -fn SP -off %f -end 0.99 -pow %f -c %f    \\\n' % (win1Val, win2Val, c))
    elif win == 'EM':
        outy.write('| nmrPipe  -fn EM -lb %f -c %f    \\\n' % (win1Val, c))
    else:
        outy.write('| nmrPipe  -fn GM -g1 %f -g2 %f -g3 0.0 -c %f    \\\n' % (win1Val, win2Val, c))

    outy.write('| nmrPipe  -fn ZF -zf 2                               \\\n')
    outy.write('| nmrPipe  -fn FT -auto                               \\\n')
    outy.write('| pipe2xyz -out %s -x -verb -ov\n' % (frame._spec_output_dir()+'/slice.ft1',))
    outy.close()

    return


def make_proc_script_2d_state(process, state: ProcessingScriptState, outfile, lp='n'):




    #Dealing with transposes is a pain.
    #ultimately we desire for manco:
    # X   Y   Z    A
    #13Ca 1Hz 13Cy HNx
    # but these start off as:
    # A    Z  Y   X
    #planes will be Z/A.
    #X and Y are the ones that will be scanned over
    #indirect dimension needs to get placed on 'A' asap.
    #The first few functions are for normal FTs.
    #by this point we get the right transposition.
    #everything is then done with ZTP and TP.
    #annoyingly, nmrpipe doens't have a 4D TP.
    #in the case of this spectrum:
    #F1180 on: Y 1Hz and A HNx

    outy=open(outfile,'w')
    outy.write('#!/bin/csh -f\n')
    outy.write('%s\n' % ('set ft4trec='+process._spec_output_dir()+'/test.fid'))

    outy.write('echo \'Processing 2D\'\n')
    outy.write('showhdr $ft4trec\n')

    outy.write('echo Processing XY Dimensions Without LP ......\n')
    outy.write('if( -e %s/test.ft2) rm -rf %s/test.ft2\n' % (process._spec_output_dir(),process._spec_output_dir()))
    outy.write('xyz2pipe -in $ft4trec -verb -x \\\n')
    outy.write('# Process X \\\n')
    if(state.checked('cb_baseSol')):
        outy.write('| nmrPipe -fn SOL \\\n')

    do_flags(process,state.control('windowBox0'),state.control('win2Val0'),state.control('win3Val0'),state.control('firstPoint0'),False,False,'n',state.control('p0'),state.control('p1'),state.control('cb_ft0'),outy)

    do_extract_f1(process, outy)


    outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
    if(lp=='y'):

        nuslist = process._current_nus_schedule()
        # SMILE controls belong to ProcessingFrame.  Script generation uses
        # the immutable ProcessingScriptState snapshot so it also works when
        # ProcessFrame no longer aliases those wx widgets (notably 2D + LP).
        try:
            maxIter = float(state.value('maxIterBox', 0))
        except (TypeError, ValueError):
            maxIter = 0.0

        dic,data=ng.pipe.read(process._spec_output_dir()+'/test.fid')
        xT=int(dic['FDF1TDSIZE'])  #35  Cy
        #aT=int(dic['FDF2TDSIZE'])  #951 Hx
        yT=int(dic['FDF3TDSIZE'])  #16  Hz
        zT=int(dic['FDF4TDSIZE'])  #35  Ca

        #AFTER ZTP
        #zT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE  #DIRECT DIMENSION
        #yT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
        #xT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
        #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE
        #print ord
        #print xT,yT,zT,aT

        #AFTER ZTP AND TP
        #zT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE  #DIRECT DIMENSION
        #xT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
        #yT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
        #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE

        #logging.info(xT,yT,zT)
        #40 0 40 512
        #### CHARLIE SET NCPUS to os.cpu_count() RATHER THAN OFF UNIDEC WINDOW
        if nuslist != '':
            raw_dir = str(process._raw_output_dir()).strip()
            outy.write('| nmrPipe  -fn SMILE  -nDim 2 -nThread %i -report 1 -sample %s/%s \\\n' % (process.ncpus,raw_dir,nuslist))
            if(maxIter>0):
                outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
            else:
                outy.write('                		       -maxIter 50  \\\n') #default is 800

            #AFTER ZIP AND TP (CHECK!!)
            do_smile(process,xT,'x',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Z
            # Higher indirect dimensions are handled by the 3D/4D builders.

        elif(state.checked('cb_lp1')): #IF NOT RECONSTRUCTING
            outy.write('| nmrPipe  -fn SMILE  -nDim 2 -nThread %i -report 1  \\\n' % (process.ncpus))
            #outy.write('                #		       -sigma 0    \\\n')
            #outy.write('                ##		       -thresh 0.1 \\\n')
            if(maxIter>0):
                outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
            else:
                outy.write('                		       -maxIter 50  \\\n') #default is 800

            do_smile(process,xT,'x',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Z

    outy.write('# Process Y \\\n')

    do_flags(process,state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),state.control('firstPoint1'),state.control('cb_f1180'),state.control('cb_lp1'),lp,state.control('p0_1'),state.control('p1_1'),state.control('cb_ft1'),outy)

    #CHARLIE DID NOT WRTAP BASE/POLY IN TRANSPOSE STATEMENTS
    if(state.checked('cb_baseLin')):
        outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
        outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
        outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
    if(state.checked('cb_basepol')):
        outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
        outy.write('| nmrPipe -fn POLY -auto\\\n')
        outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
    if(state.checked('cb_basepol1')):
        outy.write('| nmrPipe -fn POLY -auto\\\n')

    #it might be that this is needed only for bruker! makes sure 1H ends up on outermost dimension.
    outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
        
    #if(state.checked('cb_baseLin')):
    #    outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
    #if(state.checked('cb_basepol')):
    #    outy.write('| nmrPipe -fn POLY -auto\\\n')
    #if(state.checked('cb_basepol1')): #baseline along
    #    outy.write('| nmrPipe -fn ZTP \\\n')
    #    outy.write('| nmrPipe -fn POLY -auto\\\n')
    #    outy.write('| nmrPipe -fn ZTP \\\n')
    #if(state.checked('cb_basepol2')): #baseline along
    #    outy.write('| nmrPipe -fn TP \\\n')
    #    outy.write('| nmrPipe -fn POLY -auto\\\n')
    #    outy.write('| nmrPipe -fn TP \\\n')

    outy.write('| pipe2xyz -out %s -x -verb -ov\n' % (process._spec_output_dir()+'/test.ft2',))
    outy.close()

    return



def make_proc_script_2d(frame, outfile, lp='n'):
    """Legacy frame-shaped entry point retained during staged migration."""
    state = ProcessingScriptState.capture(frame, {}, PROCESSING_CONTROL_NAMES)
    return make_proc_script_2d_state(frame, state, outfile, lp=lp)

def make_proc_script_3dp_state(process, state: ProcessingScriptState, outfile, lp='n'):
    #Dealing with transposes is a pain.
    #ultimately we desire for manco:
    # X   Y   Z    A
    #13Ca 1Hz 13Cy HNx
    # but these start off as:
    # A    Z  Y   X
    #planes will be Z/A.
    #X and Y are the ones that will be scanned over
    #indirect dimension needs to get placed on 'A' asap.
    #The first few functions are for normal FTs.
    #by this point we get the right transposition.
    #everything is then done with ZTP and TP.
    #annoyingly, nmrpipe doens't have a 4D TP.
    #in the case of this spectrum:
    #F1180 on: Y 1Hz and A HNx

    outy=open(outfile,'w')
    outy.write('#!/bin/csh -f\n')
    outy.write('%s\n' % ('set ft4trec='+process._spec_output_dir()+'/fids/test%03d.fid'))

    outy.write('echo \'Processing pseudo 3D\'\n')
    outy.write('showhdr $ft4trec\n')

    outy.write('echo Processing XY Dimensions Without LP ......\n')
    outy.write('if( -e %s/test.ft2) rm -rf %s/test.ft2\n' % (process._spec_output_dir(),process._spec_output_dir()))
    outy.write('xyz2pipe -in $ft4trec -verb -x \\\n')
    outy.write('# Process X \\\n')
    if(state.checked('cb_baseSol')):
        outy.write('| nmrPipe -fn SOL \\\n')

    do_flags(process,state.control('windowBox0'),state.control('win2Val0'),state.control('win3Val0'),state.control('firstPoint0'),False,False,'n',state.control('p0'),state.control('p1'),state.control('cb_ft0'),outy)

    do_extract_f1(process, outy)

    outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
    if(lp=='y'):

        #if(state.checked('cb_lp1') and lp=='y'):
        outy.write('| nmrPipe  -fn LP -fb                 \\\n')

        """
        nuslist = process._current_nus_schedule()
        maxIter=float(state.value('maxIterBox', 0))
        try:
            maxIter=float(maxIter)
        except:
            pass

        dic,data=ng.pipe.read(process._spec_output_dir()+'/fids/test001.fid')
        xT=int(dic['FDF1TDSIZE'])  #35  Cy
        ##aT=int(dic['FDF2TDSIZE'])  #951 Hx
        yT=int(dic['FDF3TDSIZE'])  #16  Hz
        zT=int(dic['FDF4TDSIZE'])  #35  Ca

        #bbbbbbbbbbb

        #AFTER ZTP
        #zT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE  #DIRECT DIMENSION
        #yT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
        #xT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
        #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE
        #print ord
        #print xT,yT,zT,aT

        #AFTER ZTP AND TP
        #zT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE  #DIRECT DIMENSION
        #xT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
        #yT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
        #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE

        #print(xT,yT,zT)
        #sys.exit(100)
        #40 0 40 512

        if(nuslist!=''):

            outy.write('| nmrPipe  -fn SMILE  -nDim 2 -nThread %i -report 1 -sample %s/%s \\\n' % (process.ncpus,process._raw_output_dir(),nuslist))
            if(maxIter>0):
                outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
            else:
                outy.write('                		       -maxIter 50  \\\n') #default is 800

            #AFTER ZIP AND TP (CHECK!!)
            do_smile(process,xT,'x',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Z
            #do_smile(process,yT,'y',state.control('p0_2'),state.control('p1_2'),state.control('cb_f2180'),state.control('cb_lp2'),state.control('windowBox2'),state.control('firstPoint2'),state.control('win2Val2'),outy) #Y

        elif(state.checked('cb_lp1')): #IF NOT RECONSTRUCTING
            outy.write('| nmrPipe  -fn SMILE  -nDim 3 -nThread %i -report 1  \\\n' % (process.ncpus))
            #outy.write('                #		       -sigma 0    \\\n')
            #outy.write('                ##		       -thresh 0.1 \\\n')
            if(maxIter>0):
                outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
            else:
                outy.write('                		       -maxIter 50  \\\n') #default is 800

            do_smile(process,xT,'x',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Z
        """

    outy.write('# Processes Y \\\n')

    do_flags(process,state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),state.control('firstPoint1'),state.control('cb_f1180'),state.control('cb_lp1'),'n',state.control('p0_1'),state.control('p1_1'),state.control('cb_ft1'),outy) #LP

    # outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule

    if(state.checked('cb_baseLin')):
        outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
    if(state.checked('cb_basepol')):
        outy.write('| nmrPipe -fn POLY -auto\\\n')

    outy.write('| nmrPipe  -fn TP \\\n') 
        
    #outy.write('| nmrPipe  -fn ZTP \\\n')  #align this with peak liost mode
    outy.write('| pipe2xyz -out %s -x -verb -ov\n' % (process._spec_output_dir()+'/test.ft2',))

    # Projection generation is deliberately not embedded in this processing
    # script.  execute_process_script() dispatches from DatasetTopology after
    # processing and uses MakeProj3P/nmrglue for 2 spectral + 1 pseudo data.
    # Keeping the legacy projection command here would incorrectly treat test.ft2 as ordinary 2D.
    outy.close()


    return


def make_proc_script_3d_state(process, state: ProcessingScriptState, outfile, lp='n'):
    #Dealing with transposes is a pain.
    #ultimately we desire for manco:
    # X   Y   Z    A
    #13Ca 1Hz 13Cy HNx
    # but these start off as:
    # A    Z  Y   X
    #planes will be Z/A.
    #X and Y are the ones that will be scanned over
    #indirect dimension needs to get placed on 'A' asap.
    #The first few functions are for normal FTs.
    #by this point we get the right transposition.
    #everything is then done with ZTP and TP.
    #annoyingly, nmrpipe doens't have a 4D TP.
    #in the case of this spectrum:
    #F1180 on: Y 1Hz and A HNx

    outy=open(outfile,'w')

    # print process.parent.build

    outy.write('#!/bin/csh -f\n')

    # if process.parent.build == "Windows":
    #     outy.write('cd "$(dirname "$0")/../"\n')
    #     outy.write('%s\n' % ('ft4trec='+process._spec_output_dir()+'/fids/test%03d.fid'))
    # else:
    outy.write('%s\n' % ('set ft4trec='+process._spec_output_dir()+'/fids/test%03d.fid'))

    outy.write('echo \'Processing 3D\'\n')
    outy.write('showhdr $ft4trec\n')
    outy.write('ls\n')
    # if process.parent.build == "Windows":
    #     outy.write('if [ -e %s/XYZA ]; then rm -rf %s/XYZA; fi\n' %(process._spec_output_dir(),process._spec_output_dir()) )
    # else:
    outy.write('if( -e %s/XYZA ) rm -rf %s/XYZA\n' %(process._spec_output_dir(),process._spec_output_dir()) )

    outy.write('mkdir %s/XYZA\n' % process._spec_output_dir())
    outy.write('echo Processing XY Dimensions Without LP ......\n')
    # if process.parent.build == "Windows":
    #     outy.write('if [ -e %s/ft ]; then rm -rf %s/ft; fi\n' % (process._spec_output_dir(),process._spec_output_dir()))
    # else:
    outy.write('if( -e %s/ft) rm -rf %s/ft\n' % (process._spec_output_dir(),process._spec_output_dir()))
    outy.write('xyz2pipe -in $ft4trec -verb -x \\\n')
    # if process.parent.build != "Windows":
    outy.write('# Process X \\\n')

    if(state.checked('cb_baseSol')):
        outy.write('| nmrPipe -fn SOL \\\n')

    do_flags(process,state.control('windowBox0'),state.control('win2Val0'),state.control('win3Val0'),state.control('firstPoint0'),False,False,'n',state.control('p0'),state.control('p1'),state.control('cb_ft0'),outy)

    do_extract_f1(process, outy)


    #if(state.checked('cb_baseLin')):
    #    outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
    #if(process.cb_basePol.IsChecked()):
    #    outy.write('| nmrPipe -fn POLY -auto\\\n')

    if(lp=='y'):

        nuslist = process._current_nus_schedule()
        maxIter=float(state.value('maxIterBox', 0))
        try:
            maxIter=float(maxIter)
        except:
            pass

        dic,data=ng.pipe.read(process._spec_output_dir()+'/fids/test001.fid')
        xT=int(dic['FDF1TDSIZE'])  #35  Cy
        #aT=int(dic['FDF2TDSIZE'])  #951 Hx
        yT=int(dic['FDF3TDSIZE'])  #16  Hz
        zT=int(dic['FDF4TDSIZE'])  #35  Ca

        #AFTER ZTP
        #zT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE  #DIRECT DIMENSION
        #yT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
        #xT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
        #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE
        #print ord
        #print xT,yT,zT,aT

        #AFTER ZTP AND TP
        #zT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE  #DIRECT DIMENSION
        #xT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
        #yT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
        #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE

        #logging.info(xT,yT,zT)
        #40 0 40 512
        outy.write('| nmrPipe  -fn ZTP \\\n')
        if(nuslist!=''):
            #FidPath=os.path.join(process.DataStoreBox.GetValue(),process.FidPathBox.GetValue())
            raw_dir = str(process._raw_output_dir()).strip()
            outy.write('| nmrPipe  -fn TP \\\n')  #make dimensions line up with sampling schedule
            outy.write('| nmrPipe  -fn SMILE  -nDim 3 -nThread %i -report 1 -sample %s/%s \\\n' % (process.ncpus,raw_dir,nuslist))
            if(maxIter>0):
                outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
            else:
                outy.write('                		       -maxIter 50  \\\n') #default is 800

            #AFTER ZIP AND TP (CHECK!!)
            do_smile(process,xT,'x',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Z
            do_smile(process,yT,'y',state.control('p0_2'),state.control('p1_2'),state.control('cb_f2180'),state.control('cb_lp2'),state.control('windowBox2'),state.control('win2Val2'),state.control('win3Val2'),outy) #Y

            outy.write('| nmrPipe  -fn TP \\\n')

        else: #IF NOT RECONSTRUCTING
            outy.write('| nmrPipe  -fn SMILE  -nDim 3 -nThread %i -report 1  \\\n' % (process.ncpus))
            #outy.write('                #		       -sigma 0    \\\n')
            #outy.write('                ##		       -thresh 0.1 \\\n')
            if(maxIter>0):
                outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
            else:
                outy.write('                		       -maxIter 50  \\\n') #default is 800

            do_smile(process,xT,'x',state.control('p0_2'),state.control('p1_2'),state.control('cb_f2180'),state.control('cb_lp2'),state.control('windowBox2'),state.control('win2Val2'),state.control('win3Val2'),outy) #Z
            do_smile(process,yT,'y',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Y

        outy.write('| nmrPipe  -fn ZTP \\\n')

    outy.write('| nmrPipe  -fn TP -auto \\\n')
    # if process.parent.build != "Windows":
    outy.write('# Processes Y \\\n')

    do_flags(process,state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),state.control('firstPoint1'),state.control('cb_f1180'),state.control('cb_lp1'),lp,state.control('p0_1'),state.control('p1_1'),state.control('cb_ft1'),outy)

    outy.write('| pipe2xyz -out %s -x -verb -ov\n' % (process._spec_output_dir()+'/ft/t1%03d.ft3',))
    outy.write('#if( -e %s/fids ) rm -rf %s/fids\n' % (process._spec_output_dir(),process._spec_output_dir()))
    outy.write('#put direct dimension on A asap.\n')
    outy.write('showhdr %s/ft/t1001.ft3\n' % process._spec_output_dir())

    outy.write('echo Processing Z Dimension  ......\n')
    inf=process._spec_output_dir()+'/ft/t1%03d.ft3'
    ouf=process._spec_output_dir()+'/XYZA/t2.ft3'

    # if process.parent.build == "Windows":
    #     outy.write('if [ -e %s ]; then rm %s; fi\n' % (ouf,ouf))
    # else:
    outy.write('if( -e %s ) rm %s\n' % (ouf,ouf))

    outy.write('xyz2pipe -in %s -x -verb                    \\\n' % inf)
    outy.write('| nmrPipe  -fn ZTP                                \\\n')

    do_flags(process,state.control('windowBox2'),state.control('win2Val2'),state.control('win3Val2'),state.control('firstPoint2'),state.control('cb_f2180'),state.control('cb_lp2'),lp,state.control('p0_2'),state.control('p1_2'),state.control('cb_ft2'),outy)

    outy.write('| nmrPipe  -fn TP -auto \\\n')

    if(process.tp=='bruk'): #swap the 2nd and third dimension.
        outy.write('| nmrPipe  -fn TP \\\n')
        outy.write('| nmrPipe  -fn ZTP \\\n')
        outy.write('| nmrPipe  -fn TP \\\n')

    if(state.checked('cb_basepol2')): #baseline along #1st dim?
        outy.write('| nmrPipe -fn TP \\\n')
        #outy.write('| nmrPipe  -fn ZTP \\\n')
        outy.write('| nmrPipe -fn POLY -auto\\\n')
        #outy.write('| nmrPipe  -fn ZTP \\\n')
        outy.write('| nmrPipe -fn TP \\\n')
    if(state.checked('cb_basepol1')): #baseline along #2nd dim?
        outy.write('| nmrPipe -fn ZTP \\\n')
        outy.write('| nmrPipe -fn POLY -auto\\\n')
        outy.write('| nmrPipe -fn ZTP \\\n')


    if(state.checked('cb_baseLin')):  #3rd dim?
        outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
    if(state.checked('cb_basepol')): #3rd dim?
        #outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
        outy.write('| nmrPipe -fn POLY -auto\\\n')



    outy.write('> %s\n' % ouf)
    #if process.parent.build == "Windows":
    # outy.write('if( -e %s/ft) rm -rf %s/ft\n' % (process._spec_output_dir(),process._spec_output_dir()))
    # else:
    outy.write('if( -e %s/ft) rm -rf %s/ft\n' % (process._spec_output_dir(),process._spec_output_dir()))
    outy.write('echo  %s ready\n' % ouf)
    outy.write('showhdr %s\n' % ouf)

    #make_projy_3d(frame,ouf,outy)
    #make_1d_proj_3d(frame, outy)

    outy.write('mv %s %s\n' % (ouf,process._spec_output_dir()+'/test.ft3'))
    outy.write('if( -e ft) rm -rf ft\n')
    outy.write('echo Done.\n')

    outy.close()
    return

    """
    outy.write('echo Processing Y Dimension To Linear Predict .........\n')
    inf=ouf
    ouf=process._spec_output_dir()+'/XYZA/t3.ft3'
    outy.write('showhdr %s\n' % inf)
    outy.write('if( -e %s ) rm %s\n' % (ouf,ouf))
    outy.write('cat %s \\\n' % inf)
    outy.write('| nmrPipe  -fn TP                     \\\n')
    if(state.checked('cb_f1180')):
        outy.write('| nmrPipe  -fn HT -ps90-180       \\\n')
    else:
        outy.write('| nmrPipe  -fn HT                 \\\n')
    outy.write('| nmrPipe  -fn PS -inv -hdr           \\\n')
    outy.write('| nmrPipe  -fn FT -inv                \\\n')
    outy.write('| nmrPipe  -fn ZF -inv                \\\n')
    outy.write('| nmrPipe  -fn SP -inv -hdr           \\\n')

    if(state.checked('cb_lp1') and lp=='y'):
        outy.write('| nmrPipe  -fn LP -fb                 \\\n')
    if(state.checked('cb_f1180')):
        outy.write('| nmrPipe  -fn SP -off 0.5 -end 0.99 -pow 2 -c 1.0    \\\n')
        outy.write('| nmrPipe  -fn ZF -auto                               \\\n')
        if(state.checked('cb_ft1')):
            outy.write('| nmrPipe  -fn FT -neg                               \\\n')
        else:
            outy.write('| nmrPipe  -fn FT -auto                               \\\n')
        outy.write('| nmrPipe  -fn PS -p0 %s -p1 %s -di              \\\n' % (-90+float(state.value('p0_1', 0.0)),180))
    else:
        outy.write('| nmrPipe  -fn SP -off 0.5 -end 0.99 -pow 2 -c 0.5    \\\n')
        outy.write('| nmrPipe  -fn ZF -auto                               \\\n')
        outy.write('| nmrPipe  -fn FT -auto                               \\\n')
        outy.write('| nmrPipe  -fn PS -p0 %s -p1 %s -di                 \\\n' % (state.value('p0_1', 0.0) ))
    outy.write('| nmrPipe  -fn TP                     \\\n')
    outy.write('> %s\n' % ouf)
    outy.write('if( -e %s ) rm %s\n' % (inf,inf))
    outy.write('echo Done.\n')
    outy.write('echo  %s ready\n' % ouf)

    process.MakeProj3D(ouf,outy)
    """



def make_proc_script_3dp(frame, outfile, lp='n'):
    """Legacy frame-shaped pseudo-3D entry point retained during staged migration."""
    state = ProcessingScriptState.capture(frame, {}, PROCESSING_CONTROL_NAMES)
    return make_proc_script_3dp_state(frame, state, outfile, lp=lp)


def make_proc_script_3d(frame, outfile, lp='n'):
    """Legacy frame-shaped 3D entry point retained during staged migration."""
    state = ProcessingScriptState.capture(frame, {}, PROCESSING_CONTROL_NAMES)
    return make_proc_script_3d_state(frame, state, outfile, lp=lp)

def make_proc_script_4d_state(process, state: ProcessingScriptState, outfile, lp='n'):
            #Dealing with transposes is a pain.
            #ultimately we desire for manco:
            # X   Y   Z    A
            #13Ca 1Hz 13Cy HNx
            # but these start off as:
            # A    Z  Y   X
            #planes will be Z/A.
            #X and Y are the ones that will be scanned over
            #indirect dimension needs to get placed on 'A' asap.
            #The first few functions are for normal FTs.
            #by this point we get the right transposition.
            #everything is then done with ZTP and TP.
            #annoyingly, nmrpipe doens't have a 4D TP.
            #in the case of this spectrum:
            #F1180 on: Y 1Hz and A HNx

            outy=open(outfile,'w')
            outy.write('#!/bin/csh -f\n')

            if(os.path.exists(process._spec_output_dir()+'/fids/test001001.fid')==0):
                logging.info('No conversion.')
                return -1

            outy.write('%s\n' % ('set ft4trec='+process._spec_output_dir()+'/fids/test%03d%03d.fid'))

            outy.write('echo \'Processing 4D\'\n')
            outy.write('showhdr $ft4trec\n')
            outy.write('if( -e %s/XYZA ) rm -rf %s/XYZA\n' %(process._spec_output_dir(),process._spec_output_dir()) )
            outy.write('mkdir %s/XYZA\n' % process._spec_output_dir())

            outy.write('echo Processing XY Dimensions Without LP ......\n')
            outy.write('if( -e %s/ft) rm -rf %s/ft\n' % (process._spec_output_dir(),process._spec_output_dir()))
            outy.write('xyz2pipe -in $ft4trec -verb -x \\\n')
            outy.write('# Process X \\\n')
            if(state.checked('cb_baseSol')):
                outy.write('| nmrPipe -fn SOL \\\n')

            do_flags(process,state.control('windowBox0'),state.control('win2Val0'),state.control('win3Val0'),state.control('firstPoint0'),False,False,'n',state.control('p0'),state.control('p1'),state.control('cb_ft0'),outy)

            do_extract_f1(process, outy)


            #if(state.checked('cb_baseLin')):
            #    outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
            #if(state.checked('cb_basepol')):
            #    outy.write('| nmrPipe -fn POLY -auto\\\n')

            if(lp=='y'):
                nuslist = process._current_nus_schedule()
                if(not state.checked('cb_lp1') and not state.checked('cb_lp2') and not state.checked('cb_lp3') and nuslist==''):
                    logging.info('Need to have at least one dimension selected for lp')
                    return -1

                dic,data=ng.pipe.read(process._spec_output_dir()+'/fids/test001001.fid')

                ord=process.dic['FDDIMORDER'] #order of converted fid from fid.test.com
                #xT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE
                #yT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
                #zT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
                #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE
                #print ord
                #print xT,yT,zT,aT

                #FOR EIN
                #AFTER AN ATP:
                #xT=int(dic['FDF1TDSIZE'])  #35  Cy
                #zT=int(dic['FDF3TDSIZE'])  #16  Hz
                #yT=int(dic['FDF4TDSIZE'])  #35  Ca
                #aT=int(dic['FDF2TDSIZE'])  #951 Hx

                #aT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE

                maxIter=state.value('maxIterBox', 0)
                try:
                    maxIter=float(maxIter)
                except:
                    pass

                outy.write('| nmrPipe  -fn ATP \\\n')
                if(nuslist!=''):

                    if(state.checked('cb_lp1') or state.checked('cb_lp2') or state.checked('cb_lp3')):
                        logging.info('Need LP to be switched off for reconstructions')
                        logging.info('Aborting')
                        return -1

                    outy.write('| nmrPipe  -fn ZTP \\\n')
                    outy.write('| nmrPipe  -fn TP \\\n')

                    #ASSUMES NUS SAMPLING SCHEDULE COLUMNS GO:
                    #Y,Z,A when compared to fid.test.com

                    #xT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE
                    #yT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
                    #zT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
                    #aT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE
                    #print xT,yT,zT,aT
                    #After ATP,ZTP then TP:
                    #aT=dic['FDF'+str(int(ord[0]))+'TDSIZE'] #X  dim 1 FDF2TDSIZE
                    xT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
                    yT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
                    zT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE
                    #print xT,yT,zT,aT
                    #return -1

                    #FidPath=os.path.join(frame.DataStoreBox.GetValue(),frame.FidPathBox.GetValue())
                    raw_dir = str(process._raw_output_dir()).strip()
                    outy.write('| nmrPipe  -fn SMILE  -nDim 4 -nThread %i -report 1 -sample %s/%s \\\n' % (process.ncpus,raw_dir,nuslist))
                    #outy.write('                #		       -sigma 0    \\\n')
                    #outy.write('                ##		       -thresh 0.1 \\\n')
                    if(maxIter>0):
                        outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
                    else:
                        outy.write('                		       -maxIter 50  \\\n') #default is 800
                    #FAST
                    #outy.write('-maxIter 1  \\\n')
                    #outy.write('-sigma 1    \\\n')
                    #outy.write('-thresh 1 \\\n')

                    #with ATP ZTP TP:
                    do_smile(process,xT,'x',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy)#Y
                    do_smile(process,zT,'z',state.control('p0_3'),state.control('p1_3'),state.control('cb_f3180'),state.control('cb_lp3'),state.control('windowBox3'),state.control('win2Val3'),state.control('win3Val3'),outy)#A
                    do_smile(process,yT,'y',state.control('p0_2'),state.control('p1_2'),state.control('cb_f2180'),state.control('cb_lp2'),state.control('windowBox2'),state.control('win2Val2'),state.control('win3Val2'),outy)#Z

                    outy.write('| nmrPipe  -fn TP \\\n')  #reverse transposes
                    outy.write('| nmrPipe  -fn ZTP \\\n')

                else:
                    #AFTER AN ATP:
                    yT=dic['FDF'+str(int(ord[1]))+'TDSIZE'] #Y  dim 0 FDF1TDSIZE
                    zT=dic['FDF'+str(int(ord[2]))+'TDSIZE'] #Z  dim 2 FDF3TDSIZE
                    xT=dic['FDF'+str(int(ord[3]))+'TDSIZE'] #A  dim 3 FDF4TDSIZE

                    outy.write('| nmrPipe  -fn SMILE  -nDim 4 -nThread %i -report 1  \\\n' % (process.ncpus))
                    #outy.write('                #		       -sigma 0    \\\n')
                    #outy.write('                ##		       -thresh 0.1 \\\n')
                    if(maxIter>0):
                        outy.write('                		       -maxIter %i  \\\n' % (maxIter)) #default is 800
                    else:
                        outy.write('                		       -maxIter 50  \\\n') #default is 800
                    #FAST
                    #outy.write('-maxIter 1  \\\n')
                    #outy.write('-sigma 1    \\\n')
                    #outy.write('-thresh 1 \\\n')
                    #order following ATP:
                    do_smile(process,xT,'x',state.control('p0_3'),state.control('p1_3'),state.control('cb_f3180'),state.control('cb_lp3'),state.control('windowBox3'),state.control('win2Val3'),state.control('win3Val3'),outy) #A
                    do_smile(process,zT,'z',state.control('p0_2'),state.control('p1_2'),state.control('cb_f2180'),state.control('cb_lp2'),state.control('windowBox2'),state.control('win2Val2'),state.control('win3Val2'),outy) #Z
                    do_smile(process,yT,'y',state.control('p0_1'),state.control('p1_1'),state.control('cb_f1180'),state.control('cb_lp1'),state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),outy) #Y

                outy.write('| nmrPipe  -fn ATP \\\n')

            outy.write('| nmrPipe  -fn TP -auto \\\n')
            outy.write('# Processes Y \\\n')

            do_flags(process,state.control('windowBox1'),state.control('win2Val1'),state.control('win3Val1'),state.control('firstPoint1'),state.control('cb_f1180'),state.control('cb_lp1'),lp,state.control('p0_1'),state.control('p1_1'),state.control('cb_ft1'),outy)

            outy.write('| nmrPipe  -fn TP                                         \\\n')
            outy.write('| pipe2xyz -out %s -a -verb -ov\n' % (process._spec_output_dir()+'/ft/t1%03d%03d.ft4',))
            outy.write('# YZAX \\\n')
            outy.write('#if( -e %s/fids ) rm -rf %s/fids\n' % (process._spec_output_dir(),process._spec_output_dir()))
            outy.write('#put direct dimension on A asap.\n')
            outy.write('showhdr %s/ft/t1001001.ft4\n' % process._spec_output_dir())

            outy.write('echo Processing Z Dimension  ......\n')
            inf=process._spec_output_dir()+'/ft/t1%03d%03d.ft4'
            ouf=process._spec_output_dir()+'/XYZA/t2.ft4'
            outy.write('if( -e %s ) rm %s\n' % (ouf,ouf))
            outy.write('xyz2pipe -in %s -y -verb                    \\\n' % inf)
            outy.write('# ZYAX \\\n')
            do_flags(process,state.control('windowBox2'),state.control('win2Val2'),state.control('win3Val2'),state.control('firstPoint2'),state.control('cb_f2180'),state.control('cb_lp2'),lp,state.control('p0_2'),state.control('p1_2'),state.control('cb_ft2'),outy)

            outy.write('> %s\n' % ouf)
            outy.write('#ZYAX\n')
            outy.write('if( -e %s/ft) rm -rf %s/ft\n' % (process._spec_output_dir(),process._spec_output_dir()))
            outy.write('echo  %s ready\n' % ouf)
            outy.write('showhdr %s\n' % ouf)

            outy.write('echo Processing A Dimension  .........\n')
            inf=ouf
            ouf=process._spec_output_dir()+'/XYZA/t3.ft4'
            outy.write('showhdr %s\n' % inf)
            outy.write('if( -e %s ) rm %s\n' %(ouf,ouf))
            outy.write('cat %s  \\\n' % inf)
            outy.write('| nmrPipe  -fn ZTP                                \\\n')
            outy.write('# AZYX \\\n')
            do_flags(process,state.control('windowBox3'),state.control('win2Val3'),state.control('win3Val3'),state.control('firstPoint3'),state.control('cb_f3180'),state.control('cb_lp3'),lp,state.control('p0_3'),state.control('p1_3'),state.control('cb_ft3'),outy) #do Cy

            #if(state.control('cb_lp1').IsChecked() and lp=='y'):
            #    outy.write('| nmrPipe  -fn LP -fb                               \\\n')

            process.GetLabs()
            label_heads = [str(lab)[0] if lab else '' for lab in getattr(process, 'labb', [])]
            x_lab = label_heads[0] if len(label_heads) > 0 else ''
            y_lab = label_heads[1] if len(label_heads) > 1 else ''
            z_lab = label_heads[2] if len(label_heads) > 2 else ''
            a_lab = label_heads[3] if len(label_heads) > 3 else ''

            if a_lab and a_lab == x_lab:  # if x and a dim are the same (Flemming)
                # if x and a are the same nucleus... (H C C H) eg. Flemming
                outy.write('| nmrPipe  -fn ATP                     \
    ')  # ein: Hy Cy Cx Hx flem: Cy Hy Cx Hx
                outy.write('| nmrPipe  -fn TP                     \
    ')  # ein: Hy Cy Hx Cx flem: Cy Hy Hx Cx

            elif y_lab and y_lab == x_lab:  # if x and y are the same (Varian)
                #         Hx Cx Hy Cy
                # start at Hx Hy Cx Cy
                # We are Hx Hy Cy Cx
                # start at Hx Cy Hy Cx
                outy.write('| nmrPipe  -fn TP                     \
    ')  # ein: Hx Cy Cx Hy
                outy.write('| nmrPipe  -fn ATP                     \
    ')  # ein: Hy Cy Cx Hx
                outy.write('| nmrPipe  -fn TP                     \
    ')  # ein: Hy Cy Hx Cx

            elif z_lab and z_lab == x_lab:  # if x and z are the same
                # if x and a are different nuclei... (H C H C) eg. EIN
                outy.write('| nmrPipe  -fn ZTP                     \
    ')  # ein: Hx Cy Cx Hy
                outy.write('| nmrPipe  -fn ATP                     \
    ')  # ein: Hy Cy Cx Hx flem: Cy Hy Cx Hx
                outy.write('| nmrPipe  -fn TP                     \
    ')  # ein: Hy Cy Hx Cx flem: Cy Hy Hx Cx
            else:  # otherwise do nothing and sort this out manually.
                pass


            #NEEDS TESTING
            if(state.checked('cb_baseLin')):
                outy.write('| nmrPipe -fn BASE -nw 2 -nl 0% 5% 95% 100%\\\n')
            if(state.checked('cb_basepol')):
                outy.write('| nmrPipe -fn POLY -auto\\\n')
            if(state.checked('cb_basepol1')): #baseline along
                outy.write('| nmrPipe -fn ZTP \\\n')
                outy.write('| nmrPipe -fn POLY -auto\\\n')
                outy.write('| nmrPipe -fn ZTP \\\n')
            if(state.checked('cb_basepol2')): #baseline along
                outy.write('| nmrPipe -fn TP \\\n')
                outy.write('| nmrPipe -fn POLY -auto\\\n')
                outy.write('| nmrPipe -fn TP \\\n')
            #if(frame.cb_basepol3.IsChecked()): #baseline along
            #    outy.write('| nmrPipe -fn TP \\\n')
            #    outy.write('| nmrPipe -fn POLY -auto\\\n')
            #    outy.write('| nmrPipe -fn TP \\\n')



            outy.write('> %s\n' % ouf)
            outy.write('if( -e %s ) rm %s\n' % (inf,inf))
            outy.write('echo  %s ready\n' % (ouf))
            outy.write('showhdr %s\n' % (ouf))

            make_proj_4d(process,ouf,outy)
            return


def make_proc_script_4d(frame, outfile, lp='n'):
    """Legacy frame-shaped 4D entry point retained for external callers."""
    state = ProcessingScriptState.capture(frame, {}, PROCESSING_CONTROL_NAMES)
    return make_proc_script_4d_state(frame, state, outfile, lp=lp)

def write_direct_phase_script(frame, p0, p1):
    if getattr(frame, '_direct_phase_backend', lambda: 'glue')() != 'pipe':
        raise RuntimeError('Direct phase script generation is disabled for the nmrglue backend')
    outy = open(frame._direct_phase_script_path(), 'w')
    outy.write('#!/bin/csh -f\n')
    ftrec = os.path.join(frame._spec_output_dir(), 'slice.ft1')
    out = frame._direct_phased_spectrum_path()
    outy.write('set ft4trec=%s\n' % ftrec)
    outy.write("echo 'Processing 1D slice phase'\n")
    outy.write('showhdr $ft4trec\n')
    outy.write('if( -e %s) rm -rf %s\n' % (out, out))
    outy.write('nmrPipe -in $ft4trec \\n')
    outy.write('| nmrPipe -fn PS -p0 %s -p1 %s -di \\n' % (p0, p1))
    outy.write('| pipe2xyz -out %s -x -verb -ov\n' % out)
    outy.close()
    try:
        os.chmod(frame._direct_phase_script_path(), 0o755)
    except Exception:
        pass
    logging.debug(f'[process] wrote nmrPipe phase script: {frame._direct_phase_script_path()}')
    logging.debug(f'[process] phase script input={ftrec}')
    logging.debug(f'[process] phase script output={out}')
    logging.debug(f'[process] phase script p0={p0} p1={p1}')
    return frame._direct_phase_script_path()


def run_direct_phase_script(frame, p0, p1):
    if getattr(frame, '_direct_phase_backend', lambda: 'glue')() != 'pipe':
        raise RuntimeError('Direct phase script execution is disabled for the nmrglue backend')
    script = write_direct_phase_script(frame, p0, p1)
    raw_base = frame._raw_output_dir()
    spec_base = frame._spec_output_dir()
    project_root = os.path.abspath(os.path.join(spec_base, os.pardir)) if spec_base else (os.path.abspath(os.path.join(raw_base, os.pardir)) if raw_base else os.getcwd())
    abs_script = os.path.abspath(script)
    abs_input = os.path.abspath(os.path.join(spec_base, 'slice.ft1')) if spec_base else (os.path.abspath(os.path.join(raw_base, 'slice.ft1')) if raw_base else os.path.abspath('slice.ft1'))
    abs_output = os.path.abspath(frame._direct_phased_spectrum_path())

    logging.debug(f'[process] _run_direct_phase_script cwd(before)={os.getcwd()}')
    logging.debug(f'[process] _run_direct_phase_script project_root={project_root}')
    logging.debug(f'[process] _run_direct_phase_script raw_base={raw_base!r}')
    logging.debug(f'[process] _run_direct_phase_script script={script!r} abs_script={abs_script!r}')
    logging.debug(f'[process] _run_direct_phase_script input={abs_input!r} exists={os.path.exists(abs_input)}')
    logging.debug(f'[process] _run_direct_phase_script output={abs_output!r} exists={os.path.exists(abs_output)}')
    logging.debug(f'[process] _run_direct_phase_script script_exists(before)={os.path.exists(abs_script)}')
    if os.path.isdir(project_root):
        try:
            logging.debug(f'[process] _run_direct_phase_script project_root_listing={sorted(os.listdir(project_root))[:20]}')
        except Exception as exc:
            logging.debug(f'[process] _run_direct_phase_script project_root_listing FAILED: {exc!r}')

    try:
        result = subprocess.run(['csh', abs_script], cwd=project_root or None, capture_output=True, text=True, check=True)
        logging.debug(f'[process] _run_direct_phase_script cwd(after)={os.getcwd()}')
        if result.stdout:
            print(result.stdout, end='')
        if result.stderr:
            print(result.stderr, end='')
    except subprocess.CalledProcessError as exc:
        logging.debug(f'[process] _run_direct_phase_script returncode={exc.returncode}')
        if exc.stdout:
            print(exc.stdout, end='')
        if exc.stderr:
            print(exc.stderr, end='')
        logging.debug(f'[process] _run_direct_phase_script script_exists(after)={os.path.exists(abs_script)}')
        logging.debug(f'[process] _run_direct_phase_script output_exists(after)={os.path.exists(abs_output)}')
        raise RuntimeError(f'nmrPipe phase script failed with exit code {exc.returncode}') from exc
    return frame._direct_phased_spectrum_path()


def _build_script(frame, outfile: str, lp: str = 'n'):
    frame.GetSpectrometerType()
    if frame.spectral_dim_count == 1:
        return make_proc_script_1d(frame, outfile, lp=lp)
    if frame.spectral_dim_count == 2 and frame.has_pseudo_axis:
        return make_proc_script_3dp(frame, outfile, lp=lp)
    if frame.spectral_dim_count == 2:
        return make_proc_script_2d(frame, outfile, lp=lp)
    if frame.spectral_dim_count == 3:
        return make_proc_script_3d(frame, outfile, lp=lp)
    return make_proc_script_4d(frame, outfile, lp=lp)



def _sanitize_csh_script(text: str) -> str:
    """Normalize typography that is invalid or ambiguous in executable csh.

    Processing scripts are executable source, so smart quotes/dashes must never
    reach disk. This is deliberately applied at the final render/write boundary.
    """
    return (text.replace("\u201c", '"').replace("\u201d", '"')
                .replace("\u2018", "'").replace("\u2019", "'")
                .replace("\u2014", "--").replace("\u2013", "-")
                .replace("\u2026", "..."))

def _mddnmr_root() -> Path:
    """Return the bundled MDDNMR directory, expected beside the decon package."""
    return Path(__file__).resolve().parents[2] / 'mddnmr'


def _mddnmr_runtime_preamble() -> str:
    """csh code selecting the most suitable bundled MDDNMR binary at runtime.

    On Apple Silicon, prefer the testing-only binMAC_ARM_compat layout.  That
    directory contains the patched x86_64 MDDNMR executable plus its private
    x86_64 libstdc++/libgcc runtime.  A tiny shim forces execution via Rosetta
    while leaving the upstream MDDNMR wrapper unchanged.
    """
    root = shlex.quote(str(_mddnmr_root()))
    return """# MDDNMR runtime selection (MDDNMR directory is beside decon).
# Support both <SpinHub>/mddnmr/bin... and <SpinHub>/mddnmr/mddnmr/bin... .
setenv MDD_NMR {root}
if ( ! -e "$MDD_NMR/com/mddnmrParSet.sh" && -e "$MDD_NMR/mddnmr/com/mddnmrParSet.sh" ) then
    setenv MDD_NMR "$MDD_NMR/mddnmr"
endif
set _mdd_bin = ""
set _mdd_shim = ""
set _mdd_os = `uname -s`
set _mdd_arch = `uname -m`
if ( "$_mdd_os" == "Darwin" ) then
    if ( "$_mdd_arch" == "arm64" ) then
        set _mdd_compat = "$MDD_NMR/binMAC_ARM_compat"
        echo "=== MDDNMR Apple Silicon diagnostics ==="
        echo "MDD_NMR=<$MDD_NMR>"
        echo "compat_dir=<$_mdd_compat>"
        /bin/ls -ld "$MDD_NMR" "$_mdd_compat"
        /bin/ls -l "$_mdd_compat/mddnmr4pipeN" "$_mdd_compat/libstdc++.6.dylib" "$_mdd_compat/libgcc_s.1.dylib"
        echo "=== end MDDNMR diagnostics ==="
        # Test each required compatibility file independently.
        # Test each required compatibility file independently and record readiness.
        set _mdd_compat_ok = 1
        if ( ! -e "$_mdd_compat/mddnmr4pipeN" ) set _mdd_compat_ok = 0
        if ( ! -e "$_mdd_compat/libstdc++.6.dylib" ) set _mdd_compat_ok = 0
        if ( ! -e "$_mdd_compat/libgcc_s.1.dylib" ) set _mdd_compat_ok = 0
        if ( $_mdd_compat_ok == 1 ) then
            # Verify Rosetta before constructing the wrapper shim.
            arch -x86_64 /usr/bin/true >& /dev/null
            if ( $status ) then
                echo "MDDNMR: Apple Silicon compatibility runtime found, but Rosetta 2 is unavailable."
                echo "MDDNMR: install Rosetta with: softwareupdate --install-rosetta --agree-to-license"
                exit 1
            endif
            set _mdd_bin = "$_mdd_compat"
            set _mdd_shim = "/tmp/decon_mddnmr_$$"
            mkdir -p "$_mdd_shim"
            # Do not generate a #!/bin/sh wrapper from csh: csh performs history
            # expansion on ! even inside single quotes ("Event not found").
            # A symlink is sufficient here: macOS launches the x86_64 Mach-O via
            # Rosetta automatically, while mddnmrParSet.sh still finds the
            # expected command name on PATH.
            /bin/ln -sf "$_mdd_compat/mddnmr4pipeN" "$_mdd_shim/mddnmr4pipeN"
            if ( ! -x "$_mdd_shim/mddnmr4pipeN" ) then
                echo "MDDNMR: failed to create Rosetta command link at $_mdd_shim/mddnmr4pipeN"
                exit 1
            endif
            echo "MDDNMR runtime: Apple Silicon -> x86_64 MDDNMR via Rosetta compatibility runtime"
        else
            echo "MDDNMR: Apple Silicon detected, but the testing compatibility runtime is incomplete."
            echo "MDDNMR: expected these files under $_mdd_compat:"
            echo "  mddnmr4pipeN"
            echo "  libstdc++.6.dylib"
            echo "  libgcc_s.1.dylib"
            echo "MDDNMR: use the patched binMAC_ARM_compat runtime created for this test build."
            exit 1
        endif
    else if ( -x "$MDD_NMR/binMAC_86/mddnmr4pipeN" ) then
        set _mdd_bin = "$MDD_NMR/binMAC_86"
        echo "MDDNMR runtime: Intel macOS native x86_64 binary"
    endif
else if ( "$_mdd_os" == "Linux" ) then
    set _mdd_id = ""
    set _mdd_ver = ""
    if ( -r /etc/os-release ) then
        set _mdd_id = `awk -F= '/^ID=/{{gsub(/\"/,"",$2); print tolower($2)}}' /etc/os-release`
        set _mdd_ver = `awk -F= '/^VERSION_ID=/{{gsub(/\"/,"",$2); print $2}}' /etc/os-release`
    endif
    if ( "$_mdd_id" == "ubuntu" && "$_mdd_ver" =~ 24.* && -x "$MDD_NMR/binUbuntu_24_64/mddnmr4pipeN" ) then
        set _mdd_bin = "$MDD_NMR/binUbuntu_24_64"
    else if ( "$_mdd_id" == "ubuntu" && "$_mdd_ver" =~ 18.* && -x "$MDD_NMR/binUbuntu_18_64/mddnmr4pipeN" ) then
        set _mdd_bin = "$MDD_NMR/binUbuntu_18_64"
    else if ( ( "$_mdd_id" == "centos" || "$_mdd_id" == "rhel" || "$_mdd_id" == "rocky" || "$_mdd_id" == "almalinux" ) && -x "$MDD_NMR/binCentOS7/mddnmr4pipeN" ) then
        set _mdd_bin = "$MDD_NMR/binCentOS7"
    else if ( -x "$MDD_NMR/binUbuntu64Static/mddnmr4pipeN" ) then
        set _mdd_bin = "$MDD_NMR/binUbuntu64Static"
    else if ( -x "$MDD_NMR/binCentOS64Static/mddnmr4pipeN" ) then
        set _mdd_bin = "$MDD_NMR/binCentOS64Static"
    endif
else
    echo "MDDNMR: unsupported operating system $_mdd_os"
    exit 1
endif
if ( "$_mdd_bin" == "" || ! -x "$_mdd_bin/mddnmr4pipeN" ) then
    echo "MDDNMR: no usable mddnmr4pipeN binary found under $MDD_NMR for $_mdd_os/$_mdd_arch"
    exit 1
endif
if ( ! -x "$MDD_NMR/com/mddnmrParSet.sh" ) then
    echo "MDDNMR: missing or non-executable wrapper: $MDD_NMR/com/mddnmrParSet.sh"
    echo "MDDNMR: if present, restore permission with: chmod +x $MDD_NMR/com/mddnmrParSet.sh"
    exit 1
endif
# Use csh's special $path array rather than constructing PATH with colons.
# This avoids csh interpreting ':' after a variable expansion as a modifier.
if ( "$_mdd_shim" != "" ) then
    set path = ( "$_mdd_shim" "$MDD_NMR/com" "$_mdd_bin" $path )
else
    set path = ( "$_mdd_bin" "$MDD_NMR/com" $path )
endif
echo "MDDNMR binary: $_mdd_bin/mddnmr4pipeN"
""".format(root=root)


def _mdd_value(state, name, default):
    value = state.value(name, default)
    if value is None or str(value).strip() == '':
        return default
    return str(value).strip()


def _transform_smile_script_to_mddnmr(script: str, process, state: ProcessingScriptState) -> str:
    """Create an MDDNMR bridge without changing the established SMILE generator.

    For 3D NUS data MDDNMR expects the directly detected dimension to remain X
    and the two sparse dimensions to remain indirect when the intermediate pipe
    file is written.  The legacy SMILE path deliberately applies ZTP+TP before
    SMILE to align NMRPipe's in-stream SMILE dimension ordering.  That ordering
    is *not* the ordering expected by mddnmrParSet.sh.  Therefore this transform
    is intentionally MDDNMR-only: it removes that local ZTP+TP pair, consumes
    the SMILE-only header-option lines, writes the MDDNMR input, reconstructs,
    then restores the legacy SMILE reconstruction orientation (ZTP+TP) before resuming at the existing downstream TP chain. This post-adapter is 3D MDDNMR-only.

    The source SMILE script and its builders are never modified.
    """
    import re
    ndim = int(getattr(process, 'spectral_dim_count', 2) or 2)
    if getattr(process, 'has_pseudo_axis', False) and ndim == 2:
        ndim = 3
    ext = {2: 'ft2', 3: 'ft3', 4: 'ft4'}.get(ndim, 'ft2')
    spec = str(process._spec_output_dir())
    raw = str(process._raw_output_dir()).strip()
    nus = str(process._current_nus_schedule()).strip()
    if not nus:
        raise ValueError('MDDNMR reconstruction requires an NUS sampling schedule')
    sample = os.path.join(raw, nus)
    in_file = os.path.join(spec, 'mddnmr_input.' + ext)
    out_file = os.path.join(spec, 'mddnmr_reconstructed.' + ext)
    order = ' '.join(str(i) for i in range(1, max(2, ndim)))
    method = _mdd_value(state, 'mddMethodBox', 'CS').upper()
    alg = _mdd_value(state, 'mddAlgorithmBox', 'IST').upper()
    iterations = int(float(_mdd_value(state, 'mddIterBox', '100')))
    ve = 'y' if _mdd_value(state, 'mddVEBox', 'y').lower() in ('y', 'yes', 'true', '1', 'on') else 'n'
    threads = max(1, int(getattr(process, 'ncpus', 1) or 1))
    phase = ' '.join(['0'] * ndim)

    diagnostics = (
        'echo "=== MDDNMR input header ==="\n'
        + 'showhdr %s\n' % shlex.quote(in_file)
        + 'echo "=== NUS schedule diagnostics ==="\n'
        + "awk 'NF && $1 !~ /^#/ {n++; if(n==1){for(i=1;i<=NF;i++){mn[i]=$i;mx[i]=$i}} for(i=1;i<=NF;i++){if($i<mn[i])mn[i]=$i;if($i>mx[i])mx[i]=$i}} END {printf(\"points=%.0f columns=%.0f\",n,NF); for(i=1;i<=NF;i++)printf(\" col%.0f=[%g,%g]\",i,mn[i],mx[i]); printf(\"\\n\")}' " + shlex.quote(sample) + "\n"
        + 'echo "=== end MDDNMR input diagnostics ==="\n'
    )
    bridge = (
        '| nmrPipe -ov -out %s\n' % shlex.quote(in_file)
        + 'if ( $status ) exit 1\n'
        + diagnostics
        + "set _nus_points = `awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' %s`\n" % shlex.quote(sample)
        + 'echo "MDDNMR reconstruction: %s / %s, $_nus_points sampled points"\n' % (method, alg)
        + "mddnmrParSet.sh InFile=%s OutFile=%s NUS_POINTS=$_nus_points METHOD=%s MDDTHREADS=%d selection_file=%s NUS_TABLE_ORDER='%s' SRSIZE=0.08 OVLP=1 CS_alg=%s CS_niter=%d CS_VE=%s phase='%s' f180=%s\n" % (shlex.quote(in_file), shlex.quote(out_file), method, threads, shlex.quote(sample), order, alg, iterations, ve, phase, 'n' * ndim)
        + 'if ( $status ) then\n  echo "MDDNMR reconstruction failed"\n  exit 1\nendif\n'
        + 'echo "=== MDDNMR reconstructed header ==="\n'
        + 'showhdr %s\n' % shlex.quote(out_file)
        + 'nmrPipe -in %s \\\n' % shlex.quote(out_file)
        + ('| nmrPipe -fn ZTP \\\n| nmrPipe -fn TP \\\n' if ndim == 3 else '')
    )

    # SMILE command plus its continuation options (-maxIter, -xT/-xP0/-xzf,
    # -yT/..., etc.).  Consuming these options is essential: they belong to
    # nmrPipe -fn SMILE and must not become options to `nmrPipe -in`.
    smile_block = r'\| nmrPipe\s+-fn SMILE[^\n]*\\\n(?:\s+-[^\n]*\\\n)*'

    if ndim == 3:
        # In the legacy 3D SMILE script these two transposes are immediately
        # before SMILE.  They are correct for SMILE, but make MDDNMR see the
        # cropped direct dimension as one of its NUS dimensions (e.g. 53),
        # causing valid schedule coordinates such as 205 to be rejected.
        pattern = re.compile(
            r'\| nmrPipe\s+-fn ZTP\s*\\\n'
            r'\| nmrPipe\s+-fn TP\s*\\\n'
            + smile_block
        )
        transformed, count = pattern.subn(lambda m: bridge, script, count=1)
        if count != 1:
            raise RuntimeError('Could not locate the 3D SMILE ZTP/TP reconstruction block for MDDNMR')
    else:
        # Keep existing 2D/4D experimental behaviour isolated from this 3D fix.
        pattern = re.compile(smile_block)
        transformed, count = pattern.subn(lambda m: bridge, script, count=1)
        if count != 1:
            raise RuntimeError('Could not locate the SMILE reconstruction stage to replace with MDDNMR')

    lines = transformed.splitlines(True)
    lines.insert(1 if lines and lines[0].startswith('#!') else 0, _mddnmr_runtime_preamble())
    return ''.join(lines)


def render_process_script_state(process, state: ProcessingScriptState, lp: str = 'n') -> Tuple[str, str]:
    """Render from explicit process services plus immutable processing state.

    1D generation is state-native. Higher-dimensional builders continue through
    the temporary compatibility context until their staged migration.
    """
    if not isinstance(state, ProcessingScriptState):
        raise TypeError('state must be a ProcessingScriptState')
    process.GetSpectrometerType()
    if process.spectral_dim_count <= 4:
        with TemporaryDirectory(prefix='nmrpipe_script_') as tmpdir:
            preview_path = os.path.join(tmpdir, 'preview.com')
            if process.spectral_dim_count == 1:
                make_proc_script_1d_state(process, state, preview_path, lp=('y' if lp == 'm' else lp))
            elif process.spectral_dim_count == 2 and process.has_pseudo_axis:
                make_proc_script_3dp_state(process, state, preview_path, lp=('y' if lp == 'm' else lp))
            elif process.spectral_dim_count == 2:
                make_proc_script_2d_state(process, state, preview_path, lp=('y' if lp == 'm' else lp))
            elif process.spectral_dim_count == 3:
                make_proc_script_3d_state(process, state, preview_path, lp=('y' if lp == 'm' else lp))
            else:
                make_proc_script_4d_state(process, state, preview_path, lp=('y' if lp == 'm' else lp))
            text = Path(preview_path).read_text()
            if lp == 'm':
                text = _sanitize_csh_script(_transform_smile_script_to_mddnmr(text, process, state))
            else:
                text = _sanitize_csh_script(text)
            return text, pipefile_for(process)


def write_process_script_state(process, state: ProcessingScriptState, lp: str = 'n', outfile: Optional[str] = None) -> Tuple[str, str]:
    """Atomically write a script from explicit immutable processing state."""
    process.GetSpectrometerType()
    if process.spectral_dim_count <= 4:
        if outfile is None:
            outfile = script_path_for(process, lp=lp)
        target = Path(outfile)
        target.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix='nmrpipe_script_') as tmpdir:
            tmp_path = os.path.join(tmpdir, target.name)
            if process.spectral_dim_count == 1:
                make_proc_script_1d_state(process, state, tmp_path, lp=('y' if lp == 'm' else lp))
            elif process.spectral_dim_count == 2 and process.has_pseudo_axis:
                make_proc_script_3dp_state(process, state, tmp_path, lp=('y' if lp == 'm' else lp))
            elif process.spectral_dim_count == 2:
                make_proc_script_2d_state(process, state, tmp_path, lp=('y' if lp == 'm' else lp))
            elif process.spectral_dim_count == 3:
                make_proc_script_3d_state(process, state, tmp_path, lp=('y' if lp == 'm' else lp))
            else:
                make_proc_script_4d_state(process, state, tmp_path, lp=('y' if lp == 'm' else lp))
            if lp == 'm':
                Path(tmp_path).write_text(_sanitize_csh_script(_transform_smile_script_to_mddnmr(Path(tmp_path).read_text(), process, state)))
            else:
                Path(tmp_path).write_text(_sanitize_csh_script(Path(tmp_path).read_text()))
            os.replace(tmp_path, target)
        return str(target), pipefile_for(process)


def render_process_script(frame, lp: str = 'n') -> Tuple[str, str]:
    """Return the full script text and pipefile name without saving it."""
    with TemporaryDirectory(prefix='nmrpipe_script_') as tmpdir:
        preview_path = os.path.join(tmpdir, 'preview.com')
        ret = _build_script(frame, preview_path, lp=lp)
        if ret == -1:
            raise RuntimeError('Could not generate processing script')
        text = Path(preview_path).read_text()
        return text, pipefile_for(frame)


def write_process_script(frame, lp: str = 'n', outfile: Optional[str] = None) -> Tuple[str, str]:
    """Render and atomically save the script, returning its path and pipefile."""
    if outfile is None:
        outfile = script_path_for(frame, lp=lp)
    ensure_dir = getattr(frame, '_spec_output_dir', None)
    if callable(ensure_dir):
        try:
            ensure_dir()
        except Exception:
            logging.exception('Could not ensure spec output directory before writing processing script')
    text, pipefile = render_process_script(frame, lp=lp)
    out_path = Path(outfile)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + '.tmp')
    tmp_path.write_text(_sanitize_csh_script(text))
    os.replace(tmp_path, out_path)
    return str(out_path), pipefile


def execute_process_script(frame, script_path: str, lp: str = 'n', on_finish=None, title: str = 'Processing Output') -> str:
    """Compatibility bridge to the GUI execution orchestrator.

    New GUI code should import ``decon.gui.dialogs.processing.execution``
    directly.  Keeping this bridge preserves third-party callers of the
    historical processing API while allowing this module to import headlessly.
    """
    from spinDecon.gui.dialogs.processing.execution import execute_process_script as _execute
    return _execute(frame, script_path, lp=lp, on_finish=on_finish, title=title)


class nmrPipe:
    """Shared nmrPipe helper attached to the main GUI frame."""

    def __init__(self, parent=None):
        self.parent = parent

    def bind(self, parent):
        self.parent = parent
        return self

    def _frame(self, frame=None):
        if frame is not None:
            return frame
        if self.parent is None:
            raise AttributeError('nmrPipe helper is not bound to a frame')
        return self.parent

    def script_path_for(self, frame=None, lp: str = 'n') -> str:
        return script_path_for(self._frame(frame), lp=lp)

    def pipefile_for(self, frame=None) -> str:
        return pipefile_for(self._frame(frame))

    def render_process_script_state(self, process, state: ProcessingScriptState, lp: str = 'n'):
        return render_process_script_state(process, state, lp=lp)

    def RenderProcessScriptState(self, process, state: ProcessingScriptState, lp: str = 'n'):
        return self.render_process_script_state(process, state, lp=lp)

    def write_process_script_state(self, process, state: ProcessingScriptState, lp: str = 'n', outfile: Optional[str] = None):
        return write_process_script_state(process, state, lp=lp, outfile=outfile)

    def WriteProcessScriptState(self, process, state: ProcessingScriptState, lp: str = 'n', outfile: Optional[str] = None):
        return self.write_process_script_state(process, state, lp=lp, outfile=outfile)

    def render_process_script(self, frame=None, lp: str = 'n'):
        return render_process_script(self._frame(frame), lp=lp)

    def RenderProcessScript(self, frame=None, lp: str = 'n'):
        return self.render_process_script(frame=frame, lp=lp)

    def write_process_script(self, frame=None, lp: str = 'n', outfile: Optional[str] = None):
        return write_process_script(self._frame(frame), lp=lp, outfile=outfile)

    def WriteProcessScript(self, frame=None, lp: str = 'n', outfile: Optional[str] = None):
        return self.write_process_script(frame=frame, lp=lp, outfile=outfile)

    def execute_process_script(self, frame=None, script_path: str = '', lp: str = 'n', on_finish=None, title: str = 'Processing Output') -> str:
        return execute_process_script(self._frame(frame), script_path, lp=lp, on_finish=on_finish, title=title)

    def ExecuteProcessScript(self, frame=None, script_path: str = '', lp: str = 'n', on_finish=None, title: str = 'Processing Output') -> str:
        return self.execute_process_script(frame=frame, script_path=script_path, lp=lp, on_finish=on_finish, title=title)

    def run_process_script(self, frame=None, script_path: str = '', lp: str = 'n', on_finish=None, title: str = 'Processing Output') -> str:
        return self.execute_process_script(frame=frame, script_path=script_path, lp=lp, on_finish=on_finish, title=title)

    def RunProcessScript(self, frame=None, script_path: str = '', lp: str = 'n', on_finish=None, title: str = 'Processing Output') -> str:
        return self.run_process_script(frame=frame, script_path=script_path, lp=lp, on_finish=on_finish, title=title)

    def write_direct_phase_script(self, frame=None, p0=0, p1=0):
        return write_direct_phase_script(self._frame(frame), p0, p1)

    def WriteDirectPhaseScript(self, frame=None, p0=0, p1=0):
        return self.write_direct_phase_script(frame=frame, p0=p0, p1=p1)

    def run_direct_phase_script(self, frame=None, p0=0, p1=0):
        return run_direct_phase_script(self._frame(frame), p0, p1)

    def RunDirectPhaseScript(self, frame=None, p0=0, p1=0):
        return self.run_direct_phase_script(frame=frame, p0=p0, p1=p1)

    def get_xmin_xmax(self, frame=None):
        return get_xmin_xmax(self._frame(frame))

    def GetXminXmax(self, frame=None):
        return self.get_xmin_xmax(frame=frame)

    def make_proc_script_1d(self, frame, outfile, lp: str = 'n'):
        return make_proc_script_1d(frame, outfile, lp=lp)

    def make_proc_script_1d_slice(self, frame, outfile, lp: str = 'n'):
        return make_proc_script_1d_slice(frame, outfile, lp=lp)

    def make_proc_script_2d(self, frame, outfile, lp: str = 'n'):
        return make_proc_script_2d(frame, outfile, lp=lp)

    def make_proc_script_3dp(self, frame, outfile, lp: str = 'n'):
        return make_proc_script_3dp(frame, outfile, lp=lp)

    def make_proc_script_3d(self, frame, outfile, lp: str = 'n'):
        return make_proc_script_3d(frame, outfile, lp=lp)

    def make_proc_script_4d(self, frame, outfile, lp: str = 'n'):
        return make_proc_script_4d(frame, outfile, lp=lp)

    def MakeProcScript1D(self, frame, outfile, lp: str = 'n'):
        return self.make_proc_script_1d(frame, outfile, lp=lp)

    def MakeProcScript1DSlice(self, frame, outfile, lp: str = 'n'):
        return self.make_proc_script_1d_slice(frame, outfile, lp=lp)

    def MakeProcScript2D(self, frame, outfile, lp: str = 'n'):
        return self.make_proc_script_2d(frame, outfile, lp=lp)

    def MakeProcScript3Dp(self, frame, outfile, lp: str = 'n'):
        return self.make_proc_script_3dp(frame, outfile, lp=lp)

    def MakeProcScript3D(self, frame, outfile, lp: str = 'n'):
        return self.make_proc_script_3d(frame, outfile, lp=lp)

    def MakeProcScript4D(self, frame, outfile, lp: str = 'n'):
        return self.make_proc_script_4d(frame, outfile, lp=lp)
