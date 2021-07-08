#! /usr/bin/python3

import sys
import os
import json
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import argparse

OBSERVATIONS = 'Samples/observations'
FEATURES = 'Genes/features'

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', \
        required=True, \
        dest = 'input_matrix',
        help='The input matrix'
    )

    parser.add_argument('-n', '--num_clusters', \
        required=True, \
        dest = 'num_clusters',
        type=int,
        help='The number of clusters for k-means.'
    )

    parser.add_argument('-j', '--iterations', \
        required=True, \
        dest = 'iterations',
        type=int,
        help='The number of iterations for k-means.'
    )

    parser.add_argument('-d', '--dim', \
        required=True, \
        dest = 'cluster_dimension',
        choices = [OBSERVATIONS, FEATURES],
        help='Which dimension to cluster on.'
    )

    parser.add_argument('-s', '--samples', \
        required=False, \
        dest = 'samples',
        help=('A comma-delimited list of the samples to run clustering on. Without'
            ' this argument, all samples are used.'
        )
    )

    parser.add_argument('-f', '--features', \
        required=False, \
        dest = 'features',
        help=('A comma-delimited list of the genes/features to run clustering on. Without'
            ' this argument, all genes/features are used.'
        )
    )

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()

    # read the input matrix:
    working_dir = os.path.dirname(args.input_matrix)
    f = os.path.join(working_dir, args.input_matrix)
    if os.path.exists(f):
        df = pd.read_table(f, index_col=0)
    else:
        sys.stderr.write('Could not find file: %s' % f)
        sys.exit(1)

    # if a subset of samples was requested, subset the matrix:
    if args.samples:
        samples_from_mtx = set(df.columns.tolist())
        requested_sample_list = [x.strip() for x in args.samples.split(',')]
        requested_sample_set = set(requested_sample_list)
        difference_set = requested_sample_set.difference(samples_from_mtx)
        if len(difference_set) > 0:
            sys.stderr.write('Requested samples differed from those in matrix: {csv}'.format(
                csv = ','.join(difference_set)
            ))
            sys.exit(1)
        df = df[requested_sample_list]

    # if a subset of genes/features was requested, subset the matrix:
    if args.features:
        requested_feature_list = list(set([x.strip() for x in args.features.split(',')]))
        df = df.loc[requested_feature_list]
        if df.shape[0] == 0:
            sys.stderr.write('After filtering for the requested genes/features, the matrix was empty.'
                ' Was the list of genes/features specified properly?')
            sys.exit(1)

    # Drop any genes that are fully missing. Regardless of which dimension we are 
    # clustering on, we can't work with NAs. In the case of clustering on 
    # observations, this just removes one of the components of the sample-expression
    # vector. In the case of clustering on the genes/features, then we *could*
    # try to cluster using the remaining genes, but it's just simpler to remove
    # the genes outright. 
    df = df.dropna(0, how='all')

    # Fill the remaining NAs with zeros
    df = df.fillna(0)

    if df.shape[0] == 0:
        sys.stderr.write('After removing the missing values, the resulting matrix was empty.')
        sys.exit(1) 

    # now run the k-means alg.
    kmeans = KMeans(n_clusters=args.num_clusters, max_iter=args.iterations)

    # in sklearn, they expect a matrix of (n_samples, n_features). In WebMEV,
    # we take our samples/observations in columns. Hence, need to transpose if
    # the user requested clustering on the samples
    if args.cluster_dimension == OBSERVATIONS:
        df = df.T

    try:
        # the fit_transform method expects a (samples, features) orientation
        transformed = kmeans.fit(df)
    except Exception as ex:
        sys.stderr.write('Encountered an error while calculating the k-means clusters.'
            ' Exiting. The error was reported as: {ex}'.format(ex=ex))
        sys.exit(1)

    # the cluster labels (as integers). A n-length vector where the values 
    # are on {0,...,k-1} (where k is the number of clusters)
    cluster_labels = kmeans.labels_

    # the cluster centers. A matrix of shape (k,p) where p is the dimension of
    # the vectors we are clustering on (i.e. X is (n,p))
    cluster_centers = kmeans.cluster_centers_

    # run a PCA to project the clusters for visualization
    pca = PCA(n_components=2)
    X = np.vstack([df.values, cluster_centers])
    try:
        pca_projections = pca.fit_transform(X)
    except Exception as ex:
        sys.stderr.write('Encountered an error while projecting the clusters for'
            ' visualization. The error was reported as: {ex}'.format(ex=ex))
        sys.exit(1)

    n = df.shape[0]
    if (n + args.num_clusters) != X.shape[0]:
        sys.stderr.write('Encountered an error while attempting to create visualization'
            ' of clusters. Exiting.')
        sys.exit(1)

    projected_points = pca_projections[:n,:]
    projected_centroids = pca_projections[n:,:]

    # construct the output data structure:
    output = {
        'centroids': [],
        'points' : []
    }

    for i in range(projected_centroids.shape[0]):
        x,y = projected_centroids[i,:]
        output['centroids'].append(
            {
                'cluster_id': i,
                'x':x,
                'y':y
            }
        )

    for i in range(projected_points.shape[0]):
        x,y = projected_points[i,:]
        cluster_id = cluster_labels[i]
        name = df.index[i] # the sample or gene/feature name, depending on the clustering dim.
        output['points'].append(
            {
                'x':x,
                'y':y,
                'id': name,
                'cluster_id': int(cluster_id)
            }
        )
    fname = os.path.join(working_dir, 'kmeans_output.json')
    json.dump(output, open(fname, 'w'))

    op_outputs = {
        'kmeans_results': fname
    }
    json.dump(op_outputs, open(os.path.join(working_dir, 'outputs.json'), 'w'))
