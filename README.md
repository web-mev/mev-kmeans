### K-means clustering

This MEV-compatible analysis runs the K-means algorithm as implemented in scikit-learn (https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html).

Note that following the clustering, we project the data (and centroids) into two-dimensions using PCA.

The output is a JSON-format file giving the projected locations of the input data and centroids.

Note that this is intended for use in clustering expression data, but the k-means algorithm is agnostic to the application. Do note, however, that if you use this for general clustering, you need to submit your input matrix in (features, observations) format since that is our convention in the genomics context (i.e. the observation vectors correspond to the columns of an expression matrix). Hence, the data should be arranged such that columns correspond to observations and rows correspond to features.

---

### To run external of WebMeV:

Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/kmeans:v1 .`) 
- pull the docker image from the GitHub container repository (see https://github.com/web-mev/mev-kmeans/pkgs/container/mev-kmeans)

To run, enter the container in an interactive shell:
```
docker run -it -v$PWD:/work <IMAGE>
```
(here, we mount the current directory to `/work` inside the container)

Then, run the script:
```
/opt/software/run_kmeans.py \
    -i <path to expression matrix> \
    -d <"Samples/observations" | "Genes/features"> \
    -n <number of clusters as integer> \
    -j <number of iterations as integer> \
    -s <csv string of obs/samples to include> \
    -f <csv of features/genes to include>
```

Notes:
- This script expects the input matrix in tab-delimited format.
- The clustering dimension (i.e. whether you are clustering observations or features) is given by `-d` and you should use one of those two strings. 
- `-n` is the number of clusters you want to create. Obviously needs to be less than the total number of observations or features (depending on which direction you are clustering)
- `-j` is the maximum number of k-means iterations. A reasonable default is 300.
- (Optional) If you do not want to use all the observations or features, you can use `-s` and `-f` flags to subset the matrix. You can use either or both and it doesn't matter which direction you are using for clustering as this simply performs a subset of the matrix prior to running K-means. Both of these should be formatted as comma-delimited strings. 
