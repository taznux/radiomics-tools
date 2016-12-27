from sklean.decomposition import TruncatedSVD

svd = TruncatedSVD(n_components=10)
svd.fit(xtrain)
xtrain = svd.transform(xtrain)
