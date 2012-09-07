#!/usr/bin/env python
import sys, os
import pickle
import numpy as np
from msmbuilder.metrics import (LPRMSD, RMSD, Dihedral,
                                BooleanContact, AtomPairs,
                                ContinuousContact)
from msmbuilder.metrics.baseclasses import Vectorized

def add_argument(group, *args, **kwargs):
    if 'default' in kwargs:
        d = 'Default: {d}'.format(d=kwargs['default'])
        if 'help' in kwargs:
            kwargs['help'] += ' {d}'.format(d=d)
        else:
            kwargs['help'] = d
    group.add_argument(*args, **kwargs)

################################################################################

def add_basic_metric_parsers(metric_subparser):

    metric_parser_list = []

    #metrics_parsers = parser.add_subparsers(description='Available metrics to use.',dest='metric')
    rmsd = metric_subparser.add_parser('rmsd',
        description='''RMSD: Root mean square deviation over a set of user defined atoms
        (typically backbone heavy atoms or alpha carbons). To evaluate the distance
        between two structures, first they are rotated and translated with respect
        to one another to achieve maximum coincidence. This code is executed in parallel
        on multiple cores (but not multiple boxes) using OMP.''')
    add_argument(rmsd, '-a', dest='rmsd_atom_indices', help='Atom indices to use in RMSD calculation. Pass "all" to use all atoms.', 
        default='AtomIndices.dat')
    metric_parser_list.append(rmsd)

    #rmsd_subparsers = rmsd.add_subparsers()
    #rmsd_subparsers.metric = 'rmsd'

    dihedral = metric_subparser.add_parser('dihedral',
        description='''DIHEDRAL: For each frame in the simulation data, we extract the
        torsion angles for the class of angles that you request (phi/psi is recommended,
        but chi angles are available as well). Each frame is then reprented by a vector
        containing the sin and cosine of these dihedral angles. The distances between
        frames are computed by taking distances between these vectors in R^n. The
        euclidean distance is recommended, but other distance metrics are available
        (cityblock, etc). This code is executed in parallel on multiple cores (but
        not multiple boxes) using OMP. ''') 
    add_argument(dihedral, '-a', dest='dihedral_angles', default='phi/psi',
        help='which dihedrals. Choose from phi, psi, chi. To choose multiple, seperate them with a slash')
    add_argument(dihedral, '-p', dest='dihedral_p', default=2, help='p used for metric=minkowski (otherwise ignored)')
    add_argument(dihedral, '-m', dest='dihedral_metric', default='euclidean',
        help='which distance metric', choices=Dihedral.allowable_scipy_metrics)
    metric_parser_list.append(dihedral)
#    dihedral_subparsers = dihedral.add_subparsers()
#    dihedral_subparsers.metric = 'dihedral'

    lprmsd = metric_subparser.add_parser('lprmsd',
        description='''LPRMSD: RMSD with the ability to to handle permutation-invariant atoms.
        Solves the assignment problem using a linear programming solution (LP). Can handle aligning
        on some atoms and computing the RMSD on other atoms.:''')
    add_argument(lprmsd, '-a', dest='lprmsd_atom_indices', help='Regular atom indices. Pass "all" to use all atoms.', default='AtomIndices.dat')
    add_argument(lprmsd, '-l', dest='lprmsd_alt_indices', default=None,
        help='''Optional alternate atom indices for RMSD. If you want to align the trajectories
        using one set of atom indices but then compute the distance using a different
        set of indices, use this option. If supplied, the regular atom_indices will
        be used for the alignment and these indices for the distance calculation''')
    add_argument(lprmsd, '-P', dest='lprmsd_permute_atoms', default=None, help='''Atom labels to be permuted.
    Sets of indistinguishable atoms that can be permuted to minimize the RMSD. On disk this should be stored as
    a list of newline separated indices with a "--" separating the sets of indices if there are
    more than one set of indistinguishable atoms.  Use "-- (integer)" to include a subset in the RMSD (to avoid undesirable boundary effects.)''')
    metric_parser_list.append(lprmsd)
#    lprmsd_subparsers = lprmsd.add_subparsers()
#    lprmsd_subparsers.metric = 'lprmsd'

    contact = metric_subparser.add_parser('contact',
        description='''CONTACT: For each frame in the simulation data, we extract the
        contact map (presence or absense of "contacts")  between residues. Each frame is then
        represented as a boolean valued vector containing information between the presence or
        absense of certain contacts. The contact vector can either include all possible pairwise
        contacts, only the native contacts, or any other set of pairs of residues. The distance with
        which two residues must be within to classify as "in contact" is also settable, and can
        dependend on the contact (e.g. 5 angstroms from some pairs, 10 angstroms for other pairs).
        Furthermore, the sense in which the distance between two residues is computed can be
        either specified as "CA", "closest", or "closest-heavy", which will respectively compute
        ("CA") the distance between the residues' alpha carbons, ("closest"), the closest distance between any pair of
        atoms i and j such that i belongs to one residue and j to the other residue, ("closest-heavy"), 
        or the closest distance between any pair of NON-HYDROGEN atoms i and j such that i belongs to
        one residue and j to the other residue. This code is executed in parallel on multiple cores (but
        not multiple boxes) using OMP.''')
    add_argument(contact, '-c', dest='contact_which', default='all',
        help='Path to file containing 2D array of the contacts you want, or the string "all".')
    add_argument(contact, '-C', dest='contact_cutoff', default=0.5, help='Cutoff distance in nanometers. If you pass -1, then the contact "map" will be a matrix of residue-residue distances. Passing a number greater than 0 means the residue-residue distance matrix will be converted to a boolean matrix, one if the distance is less than the specified cutoff')
    add_argument(contact, '-f', dest='contact_cutoff_file', help='File containing residue specific cutoff distances (supercedes the scalar cutoff distance if present).',default=None)
    add_argument(contact, '-s', dest='contact_scheme', default='closest-heavy', help='contact scheme.',
        choices=['CA', 'closest', 'closest-heavy'])
    metric_parser_list.append(contact)
    #contact_subparsers = contact.add_subparsers()
    #contact_subparsers.metric = 'contact'
    
    atompairs = metric_subparser.add_parser('atompairs',description='''ATOMPAIRS: For each frame, we
        represent the conformation as a vector of particular atom-atom distances. Then the distance
        between frames is calculated using a specified norm on these vectors. This code is executed in
        parallel (but not multiple boxes) using OMP.''') 
    add_argument(atompairs, '-a', dest='atompairs_which',
        help='path to file with 2D array of which atompairs to use.', default='AtomPairs.dat')
    add_argument(atompairs, '-p', dest='atompairs_p', default=2, help='p used for metric=minkowski (otherwise ignored)')
    add_argument(atompairs, '-m', dest='atompairs_metric', default='cityblock',
        help='which distance metric', choices=AtomPairs.allowable_scipy_metrics)
    metric_parser_list.append(atompairs)
    #atompairs_subparsers = atompairs.add_subparsers()
    #atompairs_subparsers.metric = 'atompairs'
    
    picklemetric = metric_subparser.add_parser('custom', description="""CUSTOM: Use a custom
        distance metric. This requires defining your metric and saving it to a file using
        the pickle format, which can be done fron an interactive shell. This is an EXPERT FEATURE,
        and requires significant knowledge of the source code's architecture to pull off.""")
    add_argument(picklemetric, '-i', dest='picklemetric_input', required=True,
        help="Path to pickle file for the metric")
    metric_parser_list.append(picklemetric)
#    picklemetric_subparsers = picklemetric.add_subparsers()
#    picklemetric_subparsers.metric = 'custom'
    
    return metric_parser_list

    ################################################################################

def add_layer_metric_parsers(metric_subparser):

    #layer_metrics = parser.add_subparsers( description='''Available Hierarchical Metrics to use. Note, 
    #    you must also specify a basic metric to prepare the trajectory with. For example if you used
    #    tICA on dihedrals you would do something like "tica --pca PCAObject.h5 --nv 10 dihedral -a phi/psi"''')

    tica = metric_subparser.add_parser( 'tica', description='''
        TICA: This metric is based on a variation of PCA which looks for the slowest d.o.f.
        in the simulation data. See (Schwantes, C.R., Pande, V.S. In Prep.) for details. Or
        contact Christian at schwancr@stanford.edu for an explanation. In addition to these 
        You must provide an additional metric you used to prepare the trajectories in the
        training step.''')

    required = tica.add_argument_group('required')
    choose_one = tica.add_argument_group('selecting projection vectors (choose_one)')

    add_argument(tica,'-p',dest='p',help='p value for p-norm')
    add_argument(tica,'-m',dest='projected_metric',help='metric to use in the projected space',
        choices= Vectorized.allowable_scipy_metrics )
    add_argument(required, '--po','--projection', dest='proj_object', help='tICA Object which was prepared by tICA_train.py')
    add_argument(choose_one, '--nv', dest='num_vecs', help='Choose the top <-n> eigenvectors based on their eigenvalues')
    add_argument(choose_one, '--ab',dest='abs_min', help='Choose all eigenvectors with eigenvalues grater than <--ab>.') 
    add_argument(choose_one, '--ev',dest='expl_var', help='Choose eigenvectors so that their eigenvalues account for <--ev> percent of the total "variance". Note that this really only makes sense when doing PCA, where the total variance is the sum of the eigenvalues.') 
    tica.metric_parser_list = []
    tica_subparsers = tica.add_subparsers( dest='sub_metric', description='''  
        Available metrics to use in preparing the trajectory before projecting.''' )
    tica.metric_parser_list = add_basic_metric_parsers(tica_subparsers)

    return tica.metric_parser_list

def add_metric_parsers(parser, add_layer_metrics=False):

    metric_parser_list = []

    metric_subparser = parser.add_subparsers( dest='metric', description ='Available metrics to use.' )

    if add_layer_metrics:
        metric_parser_list.extend( add_layer_metric_parsers(metric_subparser) )

    metric_parser_list.extend( add_basic_metric_parsers(metric_subparser) )

    return metric_parser_list

def construct_basic_metric(metric_name, args):
    if metric_name == 'rmsd':
        if args.rmsd_atom_indices != 'all':
            atom_indices = np.loadtxt(args.rmsd_atom_indices, np.int)
        else:
            atom_indices = None
        metric = RMSD(atom_indices)#, omp_parallel=args.rmsd_omp_parallel)

    elif args.metric == 'dihedral':
        metric = Dihedral(metric=args.dihedral_metric,
            p=args.dihedral_p, angles=args.dihedral_angles)
             
    elif args.metric == 'lprmsd':
        if args.lprmsd_atom_indices != 'all':
            atom_inds = np.loadtxt(args.lprmsd_atom_indices, dtype=np.int)
        else:
            atom_inds = None

        if args.lprmsd_permute_atoms is not None:
            permute_inds = ReadPermFile(args.lprmsd_permute_atoms)
        else:
            permute_inds = None

        if args.lprmsd_alt_indices is not None:
            alt_inds = np.loadtxt(args.lprmsd_alt_indices, np.int)
        else:
            alt_inds = None
    
        metric = LPRMSD(atom_inds, permute_inds, alt_inds)
         
    elif metric_name == 'contact':
        if args.contact_which != 'all':
            contact_which = np.loadtxt(args.contact_which,np.int)
        else:
            contact_which = 'all'

        if args.contact_cutoff_file != None: #getattr(args, 'contact_cutoff_file'):
            contact_cutoff = np.loadtxt(args.contact_cutoff_file, np.float)            
        elif args.contact_cutoff != None:
            contact_cutoff = float( args.contact_cutoff )
        else:
            contact_cutoff = None
             
        if contact_cutoff != None and contact_cutoff < 0:
            metric = ContinuousContact(contacts=contact_which,
                scheme=args.contact_scheme)
        else:
            metric = BooleanContact(contacts=contact_which,
                cutoff=contact_cutoff, scheme=args.contact_scheme)
     
    elif metric_name == 'atompairs':
        if args.atompairs_which != None:
            pairs = np.loadtxt(args.atompairs_which, np.int)
        else:
            pairs = None

        metric = AtomPairs(metric=args.atompairs_metric, p=args.atompairs_p,
            atom_pairs=pairs)
             
    elif metric_name == 'custom':
        with open(args.picklemetric_input) as f:
            metric = pickle.load(f)
            print '#'*80
            print 'Loaded custom metric:'
            print metric
            print '#'*80
    else:
        raise Exception("Bad metric")
     
    return metric

def construct_layer_metric(metric_name, args ):
    if metric_name == 'tica':
        sub_metric = construct_basic_metric( args.sub_metric, args )

        return metrics_projection.RedDimPNorm( args.proj_object, prep_with = sub_metric, num_vecs = args.num_vecs, abs_min = args.abs_min, metric = args.metric, p = args.p )

def construct_metric( args ):
    if hasattr( args, 'sub_metric' ):
        return construct_layer_metric( args.metric, args )
    else:
        return construct_basic_metric( args.metric, args )