"rm.matrix.validation" <-
  function(rand.mat, unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51, discard.outliers = TRUE)
  {
    if(!is.matrix(rand.mat))  stop("\n\n rm.matrix.validation: 'rand.mat' must be a matrix.\n\n")
    if(!is.numeric(rand.mat)) stop("\n\n rm.matrix.validation: 'rand.mat' must be a numeric matrix.\n\n")
    if(!unfold.method %in% c("gaussian", "spline")) stop("\n\n  rm.matrix.validation: parameter 'unfold.method' must be one of 'gaussian' or 'spline'.\n\n")

    N = nrow(rand.mat)
    M = ncol(rand.mat)
    cat(paste("\n ", N, "times", M, "matrix read.\n\n"))
    cat(paste("  Unfolding:", unfold.method, "\n"))
    if (unfold.method == "gaussian") cat(paste("  Bandwidth:", bandwidth, "\n"))
    if (unfold.method == "spline")   cat(paste("  Number of fit points:", nr.fit.points, "\n"))

    if(nrow(rand.mat) < 100) cat("\n   WARNING: 'rand.mat' seems to be quite small for this approach.\n\n")
    if(N != M) stop("\n rm.matrix.validation: 'rand.mat' must be quadratic.\n\n")
    if(!isSymmetric(rand.mat)) stop("\n rm.matrix.validation: 'rand.mat' must be symmetric.\n\n")
    if(is.complex(rand.mat)) stop("\n rm.matrix.validation: 'rand.mat' must be real-valued.\n\n")

    nr.nonzero = sum(rand.mat != 0)
    perc.nonzero = nr.nonzero/N/N*100
    sparseness = rm.get.sparseness(rand.mat)
    rank.mat = Matrix::rankMatrix(rand.mat)[1]

    cat(paste("  Number of non-zero matrix elements:", nr.nonzero, "(", perc.nonzero, "% )\n"))
    cat(paste("  Sparseness:", sparseness, "\n"))
    cat(paste("  Rank:", rank.mat, "\n"))

    fn.uf = "fit.unfold.png"
    res <- rm.ev.unfold(rand.mat, unfold.method = unfold.method, bandwidth = bandwidth, nr.fit.points = nr.fit.points,
                        discard.outliers = discard.outliers, fn = fn.uf, pop.up = FALSE, silent = TRUE)

    cat("\n")
    eigenvalues = res$eigenvalues
    ev.spacing = res$ev.spacing
    ## Results:
    results = list()
    results[["sparseness"]] = sparseness
    results[["rank"]] = rank.mat
    results[["nr.outliers.removed"]] = res$nr.outliers.removed
    return(results)
  }


## Main function: loop through thresholds and record distance of NNSD to both limiting distributions (Wigner-Dyson & Exponential)
"rm.get.threshold" <-
  function(rand.mat, nr.thresholds = 51, unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51,
           dist.method = "LL", nr.breaks = 51, discard.outliers = TRUE, discard.zeros = FALSE,
           min.mat.dim = 40, max.ev.spacing = 3, interval = NULL, interactive = TRUE, smooth.par = 0.5,
           plot.comp = TRUE, save.fit = FALSE, plot.spacing = FALSE, wait.seconds = 0)
  {
    if(!unfold.method %in% c("gaussian", "spline")) stop("\n\n  rm.get.threshold: parameter 'unfold.method' must be 'gaussian' or 'spline'.\n\n")
    if(!dist.method %in% c("LL", "KLD")) stop("\n\n  rm.get.threshold: parameter 'dist' must be 'LL' (log likelihood) or 'KLD' (Kullback-Leibler).\n\n")

    min.cell = min(abs(rand.mat[upper.tri(rand.mat, diag = F)]))
    max.cell = max(abs(rand.mat[upper.tri(rand.mat, diag = F)]))

    if(!is.null(interval)) {
      if(!is.numeric(interval) | (length(interval) != 2)) stop("\n\n rm.get.threshold: 'interval' must be a two-component numeric vector.\n\n")
      request  = ((min(interval) >= min.cell) || (min(interval) == 0)) && (max(interval) <= max.cell)
      if(!request) {
        cat(paste("\n  Range of the absolute values of the matrix elements:", signif(min.cell,5), " to ", signif(max.cell,5), "\n\n"))
        stop("\n  rm.get.threshold: parameter 'interval' must be inside the range of the absolute values of the matrix elements.\n\n")
      }
      thresholds = seq(min(interval), max(interval), len = nr.thresholds)
    } else {
      thresholds = seq(min.cell, max.cell, len = nr.thresholds)
    }

    N = nrow(rand.mat)
    cat(paste("\n ", N, "times", N, "symmetric matrix read.\n\n"))

    results = list()
    results[["unfold.method"]] = unfold.method
    results[["dist.method"]] = dist.method
    results[["tested.thresholds"]] = numeric(0)
    results[["dist.Wigner"]] = numeric(0)
    results[["dist.Expon"]] = numeric(0)
    results[["nr.zeros"]] = integer(0)
    results[["nr.uniq.ev"]]	= integer(0)
    results[["max.ev.mult"]] = integer(0)
    results[["nr.spacings"]] = integer(0)
    results[["nr.small.spacings"]] = integer(0)
    results[["perc.small.spacings"]] = integer(0)
    results[["eff.dimension"]] = integer(0)
    if(plot.comp) results[["comparison.plots"]] = character(0)
    if(save.fit) results[["cumfit.plots"]]      = character(0)
    if(plot.spacing) results[["space.plots"]]   = character(0)
    if(discard.zeros) results[["rm.dimension"]] = integer(0)
    if(discard.outliers) results[["nr.outliers.removed"]] = integer(0)
    results[["p.ks"]] = integer(0)
    results[["sse.exp"]] = integer(0)


    for (i in 1:nr.thresholds) {
      thres = thresholds[i]
      diagon = diag(rand.mat)
      rand.mat[which(abs(rand.mat) < abs(thres), arr.ind = T)] = 0
      diag(rand.mat) = diagon

      eff.mat = rm.discard.zeros(rand.mat, silent = T)

      if(discard.zeros) rand.mat = eff.mat

      if(save.fit) fn.fit = paste("RMT.Fit", i, "png", sep=".") else fn.fit <- NULL
      res <- rm.ev.unfold(rand.mat, unfold.method = unfold.method, bandwidth = bandwidth, nr.fit.points = nr.fit.points, discard.outliers = discard.outliers, fn = fn.fit, pop.up = FALSE, silent = TRUE)
      ev.spacing = res$ev.spacing


      ## Kolmogorov-Smirnov test:
      # p.val.ks.test = suppressWarnings(ks.test(ev.spacing, 'pexp', 1)$p.value)	# warnings might be generated because of ties
      p.val.ks.test = ks.test(unique(ev.spacing), 'pexp', 1)$p.value				# unique(ev.spacing) removes ties
      results[["tested.thresholds"]][i] = thres
      results[["p.ks"]][i] = p.val.ks.test


    }  	# End of loop


    thresholds   = results[["tested.thresholds"]] 		# abbreviations

    return(results)
  }



## Unfold eigenvalue spacing distribution using spline function fit
"rm.ev.unfold" <-
  function(rand.mat, unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51, discard.outliers = TRUE, fn = NULL, pop.up = FALSE, silent = TRUE)
  {

    "remove.outliers" <- function(x, factor = 1.5) {
      q25 = quantile(x, probs = 0.25); q75 = quantile(x, probs = 0.75)
      iqr = unname(q75 - q25)
      lower.threshold = q25 - (iqr * factor); upper.threshold = q75 + (iqr * factor)
      return(x[(x >= lower.threshold) & (x <= upper.threshold)])
    }

    eigenvalues = eigen(rand.mat, only.values = T)$values
    eigenvalues = eigenvalues[order(eigenvalues)]/max(abs(eigenvalues))

    ## Remove eigenvalue outliers:
    if(discard.outliers) {
      orig.number.ev = length(eigenvalues)
      eigenvalues = remove.outliers(unique(eigenvalues))
      new.number.ev = length(eigenvalues)
      nr.outliers.removed = orig.number.ev - new.number.ev
    } else {
      nr.outliers.removed = 0
    }
    if(!silent) cat(paste("  Number of discarded outlier eigenvalues:", nr.outliers.removed, "\n"))


    if(unfold.method == "gaussian") {	## Use density function obtained by Gaussian broadening
      uf <- rm.unfold.gauss(eigenvalues, bandwidth = bandwidth, fn = fn, pop.up = pop.up, silent = TRUE)
    }

    if(unfold.method == "spline") {		## Use cumulative distribution function together with cubic spline interpolation
      nr.fit.points = min(nr.fit.points, floor(0.75*length(eigenvalues)))	## experimental: the fit must be smooth, must NOT be too close to the step function (cumulative density)
      uf <- rm.unfold.spline(eigenvalues, nr.fit.points = nr.fit.points, fn = fn, pop.up = pop.up)
    }

    ## Results:
    results = list()
    results[["unfold.method"]] = unfold.method
    results[["eigenvalues"]] = uf$eigenvalues
    if(unfold.method == "gaussian") results[["bandwidth"]]  = uf$bandwidth
    if(unfold.method == "spline") results[["nr.fit.points"]] = uf$nr.fit.points
    results[["unfolded.ev"]] = uf$unfolded.ev
    results[["ev.spacing"]]  = uf$ev.spacing
    results[["plot"]] = uf$plot
    results[["nr.outliers.removed"]] = nr.outliers.removed
    return(results)
  }





"rm.unfold.gauss" <-
  function(eigenvalues, bandwidth = "nrd0", fn = NULL, pop.up = FALSE, silent = TRUE)
  {
    if(!is.numeric(eigenvalues)) stop("\n\n  rm.unfold.density: 'eigenvalues' must be a numeric vector.\n\n")
    if(!is.null(fn)) if(rm.get.file.extension(fn) != "png") stop("\n\n  rm.unfold.density: filename must have extension '.png'. \n\n")

    dens = density(eigenvalues, bw = "nrd0", adjust = 1, n = 512, kernel = "gaussian", give.Rkern = FALSE)
    bandw = dens$bw
    if(!silent) cat(paste("  Bandwith for Gaussian broadening =", signif(bandw,5), "\n"))

    midpoints <- function(x) return(x[-length(x)] + 0.5*diff(x))
    scale.function = approx(dens$x, dens$y, xout = midpoints(eigenvalues))

    ev.spacing = diff(eigenvalues)
    ev.spacing = ev.spacing * scale.function$y
    ev.spacing = ev.spacing / mean(ev.spacing)

    uf.ev = cumsum(ev.spacing)
    uf.ev = uf.ev/max(uf.ev)

    results = list()
    results[["eigenvalues"]] = eigenvalues
    results[["bandwidth"]] = bandw
    results[["unfolded.ev"]] = uf.ev
    results[["ev.spacing"]]  = ev.spacing

    return(results)
  }









"rm.get.distance" <-
  function(ev.spacing, dist.method = "LL", nr.breaks = 51)
  {
    if(!dist.method %in% c("LL", "KLD")) stop("\n\n  rm.get.distance: parameter 'dist.method' must be 'LL' (log likelihood) or 'KLD' (Kullback-Leibler).\n\n")

    results = list()
    results[["dist.method"]] = dist.method
    if(dist.method == "KLD") results[["nr.breaks"]] = nr.breaks

    if(dist.method == "LL") {						## Log Likelihood (per eigenvalue)
      evs = ev.spacing[ev.spacing != 0]
      N = length(evs)
      log.LE = -sum(evs)/N
      log.LW = log(pi/2) + sum(log(evs))/N - 0.25*pi*sum(evs^2)/N
      results[["dist.Wigner"]] = log.LW 			## log likelihood when Wigner distribution is assumed
      results[["dist.Expon"]]  = log.LE			## log likelihood when Exponential distribution is assumed
    }

    if(dist.method == "KLD") {						## Kullback-Leibler Divergence
      hsp <- hist(ev.spacing, breaks = seq(min(ev.spacing), max(ev.spacing), len = nr.breaks), plot = F)
      res.dist = kb.distance(hsp)
      results[["dist.Wigner"]] = res.dist$klw		## Kullback-Leibler divergence  observed -- Wigner surmise
      results[["dist.Expon"]]  = res.dist$klp		## Kullback-Leibler divergence  observed -- Exponential distribution
    }

    return(results)
  }





"rm.sse" <-
  function(ev.spacing, bandwidth = "nrd0", nr.points = 1000, N = 20)  # subdivide into N section with equal area
  {
    dens = density(ev.spacing, bw = bandwidth, kernel = "gaussian", n = 512)
    x = seq(min(ev.spacing), max(ev.spacing), len = nr.points)
    observed = approx(dens$x, dens$y, xout = x)$y					# calculate density at the support points
    expected = exp(-x)
    A = exp(-min(ev.spacing)) - exp(-max(ev.spacing))

    xs <- numeric(N+1)
    xs[1] = min(ev.spacing)
    for (i in 1:N) xs[i+1] = -log(exp(-xs[i]) - A/N)

    area = numeric(N) 	## area under density curve (for observed) for each section
    for (i in 1:N) {
      xsec = x[(x > xs[i]) & (x < xs[i+1])]
      xsec = c(xs[i], xsec, xs[i+1])
      ysec = approx(dens$x, dens$y, xout = xsec)$y
      area[i] = rm.trapez.int(xsec, ysec)
    }

    SSE = sum((area[i] - A/N)^2)
    return(SSE)
  }














"rm.get.sparseness" <-
  function(mat)
  {
    if(!is.matrix(mat)) stop("\n\n  rm.get.sparseness: argument 'mat' must be a matrix.\n\n")
    nr.zeros = sum(mat == 0)
    nr.cells = nrow(mat)*ncol(mat)
    sparseness = signif(nr.zeros/nr.cells,4)
    return(sparseness)
  }



"rm.trapez.int" <-
  function(x, y)
  {
    if(length(x) != length(y)) stop("\n\n  rm.trapez.int: vectors 'x' and 'y' must have the same length.\n\n")
    ind = 2:length(x)
    integral = as.double((x[ind] - x[ind-1]) %*% (y[ind] + y[ind-1]))/2
    return(integral)
  }


"wigner.surmise" <- function(x) pi/2*x*exp(-pi/4*x^2)

"rm.exp.distrib" <- function(x) exp(-x)

"wigner.semi.circle" <- function(x) 2/pi*sqrt(1-x^2)




"kb.distance" <-
  function(histo)
  {
    if(class(histo) != "histogram") stop("\n\n  kb.distance: 'histo' must be output of 'hist' function.\n\n")
    observed = histo$density
    expected.Wigner  = wigner.surmise(histo$mids)
    expected.Poisson = rm.exp.distrib(histo$mids)

    klw = kld(observed, expected.Wigner,  plot = NULL)$dist   # Kullback-Leibler divergence  observed -- Wigner surmise
    klp = kld(observed, expected.Poisson, plot = NULL)$dist   # Kullback-Leibler divergence  observed -- Exponential distrib.
    return(list(klw = klw, klp = klp))
  }





## Kullback-Leibler divergence
"kld" <-
  function(observed, expected, plot = NULL)
  {
    if(!is.numeric(expected)) stop("\n\n kld: Vector 'expected' must be numeric.\n\n")
    if(!is.numeric(observed)) stop("\n\n kld: Vector 'observed' must be numeric.\n\n")
    if(all(expected == 0)) stop("\n\n kld: All expected values are zero.\n\n")
    if(all(observed == 0)) stop("\n\n kld: All observed values are zero.\n\n")
    if(any(expected < 0)) stop("\n\n kld: expected frequency below zero.\n\n")
    if(any(observed < 0)) stop("\n\n kld: observed frequency below zero.\n\n")

    result = list()
    result[["observed"]] = observed
    result[["expected"]] = expected

    if(!is.null(plot)) {
      mtxt = "Distribution of observed and expected values"
      matplot(cbind(observed, expected), col = c("red", "blue"), main = mtxt, font.main = 1, type = "b", pch = c(1,17), ylab = "PDF", xlab = "index")
      legend("topright", c("observed", "expected"), col = c("red", "blue"), lwd = 2, pch = c(1,17))
      dev.copy(png, file = plot); dev.off()
      result[["plot"]] = plot
    }

    ind = which(observed <= 0)	# The Kullback-Leibler divergence is defined only if observed == 0 implies expected == 0
    if(length(ind) > 0) {
      observed = observed[-ind]
      expected = expected[-ind]
    }

    expected = expected/sum(expected)  	# normalize
    observed = observed/sum(observed)

    distance = sum(observed * log(observed/expected))	# Inf if any(expected == 0)

    result[["dist"]] = distance
    return(result)
  }





## Apply threshold to a matrix, the diagonal is not touched if keep.diag = T
"rm.denoise.mat" <-
  function(mat, threshold, keep.diag = TRUE)
  {
    if(!is.matrix(mat)) stop("\n  rm.denoise.mat: argument 'mat' must be a matrix\n\n")
    nr.nonzeros.before = sum(mat != 0)
    if(keep.diag) diagon = diag(mat)
    mat[which(abs(mat) < abs(threshold), arr.ind=T)] = 0
    if(keep.diag)  diag(mat) = diagon
    nr.nonzeros.after = sum(mat != 0)
    cat(paste("\n  Number of non-zero matrix elements reduced from", nr.nonzeros.before, "to", nr.nonzeros.after, "\n\n"))
    return(mat)
  }



## Remove rows and columns consisting of zeros only
"rm.discard.zeros" <-
  function(mat, tol = 0, silent = FALSE)
  {
    if(!is.matrix(mat)) stop("\n\n  discard.zero.rows: argument 'mat' must be a natrix.\n\n")
    is.null.vector <- function(x) ifelse(all(abs(x) <= tol), TRUE, FALSE)
    diagon = diag(mat)
    diag(mat) = 0
    zero.rows = apply(mat, 1, is.null.vector)
    zero.cols = apply(mat, 2, is.null.vector)
    mat = mat[!zero.rows, !zero.cols]
    diagon = diagon[!zero.rows]
    diag(mat) = diagon
    if(!silent) if((sum(zero.rows) > 0)|(sum(zero.cols) > 0)) cat(paste("  ", sum(zero.rows), "zero-rows and", sum(zero.cols), "zero-columns removed.\n"))
    return(mat)
  }



"create.rand.mat" <-
  function(size = 1000, distrib = c("norm", "unif"), mean = 0, stddev = 1)
  {
    distrib <- match.arg(distrib, c("norm", "unif"))
    cat(paste("\n  Required distribution of matrix elements:", distrib, "\n"))

    if((distrib != "norm") & ((mean!=0) | (stddev!=1))) cat("\n  Parameters 'mean' and 'stddev' ignored, apply only for the Normal distribution.\n\n")
    if(distrib == "norm") data = rnorm(size*size, mean = mean, sd = stddev)
    if(distrib == "unif") data = runif(size*size, min = -1, max = 1)

    rand.mat = matrix(data, nrow=size)
    if(distrib == "norm") rand.mat = (rand.mat + t(rand.mat))/sqrt(2)	# make symmetric, normalize
    if(distrib == "unif") rand.mat = (rand.mat + t(rand.mat))/2		# make symmetric, normalize

    m.diag  = mean(diag(rand.mat))
    sd.diag = sd(diag(rand.mat))
    ut = rand.mat[upper.tri(rand.mat)]
    m.ut = mean(ut)
    sd.ut = sd(ut)

    cat(paste("  The mean of the main diagonal is", signif(m.diag, 4), "\n"))
    cat(paste("  The std. deviation of the main diagonal is", signif(sd.diag, 4), "\n"))
    cat(paste("  The mean of the upper triangle is", signif(m.ut, 4), "\n"))
    cat(paste("  The std. deviation of the upper triangle is", signif(sd.ut, 4), "\n"))
    cat("  The matrix is real and symmetric.\n\n")

    results = list()
    results[["mean.diag"]] = signif(m.diag, 4)
    results[["stddev.diag"]] = signif(sd.diag, 4)
    results[["mean.triangle"]] = signif(m.ut, 4)
    results[["stddev.triangle"]] = signif(sd.ut, 4)
    results[["rand.matrix"]] = rand.mat
    return(results)
  }










"add.Gaussian.noise" <-
  function(mat, mean = 0, stddev = 1, symm = TRUE)
  {
    if(!is.matrix(mat)) stop("\n\n  add.Gaussian.noise: argument 'mat' must be a matrix.\n\n")
    cat(paste("\n ", nrow(mat), "x", ncol(mat), "matrix read.\n"))
    if((symm == TRUE) && (nrow(mat) != ncol(mat))) stop("\n\n  add.Gaussian.noise: matrix not quadratic - cannot be symmetrized.\n\n")

    nr.nonzero = sum(mat != 0)
    cat(paste("  This matrix contains", nr.nonzero, "non-zero cells.\n\n"))

    noise = matrix(rnorm(nrow(mat)*ncol(mat), mean = mean, sd = stddev), nrow=nrow(mat))
    mat = mat + noise
    if(symm) mat = (mat + t(mat))/sqrt(2)

    cat(paste("  Gaussian noise with mean", mean, "and standard deviation", stddev, "added.\n"))
    cat(paste("  The noise ranges from", signif(min(noise),4), "to", signif(max(noise),4), "\n"))
    if(symm) cat("  The output matrix has been symmetrized.\n\n") else cat("\n")
    return(mat)
  }






## Reorder eigenvalues and eigenvectors accordingly
"rm.reorder.ev" <-
  function(eigenvalues, eigenvec = NULL)
  {
    if(!is.null(eigenvec)) {
      if(!is.matrix(eigenvec)) stop("\n  rm.reorder.ev: 'eigenvec' must be a matrix (eigenvectors in columns)\n\n")
      if(ncol(eigenvec) != length(eigenvalues))  stop("\n  rm.reorder.ev: number of cols of 'eigenvec' must be equal to length of 'eigenvalues'\n\n")
    }
    new.order = order(eigenvalues)		# smallest first
    if(!is.null(eigenvec)) eigenvec = eigenvec[, new.order]
    eigenvalues = eigenvalues[new.order]
    results = list()
    results[["eigenvalues"]] = eigenvalues
    if(!is.null(eigenvec)) results[["eigenvec"]] = eigenvec
    return(results)
  }




"rm.connections" <-
  function(mat, nr.list = 30, abs.val = TRUE, fn = NULL)
  {
    if(!is.matrix(mat)) stop("\n\n  rm.connections: 'mat' must be a matrix.\n\n")
    if(isSymmetric(mat)) cat("\n  Matrix is symmetric.\n\n") else cat("\n  Matrix is not symmetric.\n\n")
    if(isSymmetric(mat))
      if(nr.list > sum(upper.tri(mat))) stop("\n\n  rm.connections: desired number of list elements ('nr.list') exceeds number of unique matrix elements.\n\n")
    if(!isSymmetric(mat))
      if(nr.list > nrow(mat)*ncol(mat)) stop("\n\n  rm.connections: desired number of list elements ('nr.list') exceeds number of matrix elements.\n\n")
    if(!is.null(fn))
      if(rm.get.file.extension(fn) != "txt") stop("\n\n  rm.connections: output file must have extension 'txt' (can only be saved as text file).\n\n")

    if(abs.val) temp.mat = abs(mat) else temp.mat = mat
    if(isSymmetric(temp.mat)) {
      v = temp.mat[upper.tri(temp.mat)]
      value = mat[upper.tri(mat)]
      rows = integer(0); for (i in 1:(nrow(temp.mat)-1)) rows = c(rows, 1:i)
      cols = integer(0); for (i in 2:ncol(temp.mat)) cols = c(cols, rep(i,i-1))
    } else {
      v = as.vector(temp.mat)
      value = as.vector(mat)
      rows = rep(1:nrow(temp.mat), ncol(temp.mat))
      cols = rep(1:ncol(temp.mat), each = nrow(temp.mat))
    }
    df = data.frame(row = rows, col = cols, v = v, value = value)
    df = df[order(df$v, decreasing = T),]
    df$v <- NULL
    df = df[1:nr.list,]
    rownames(df) = 1:nrow(df)
    if(!is.null(rownames(mat))) df = cbind(df, row.name = rownames(mat)[df$row])
    if(!is.null(colnames(mat))) df = cbind(df, col.name = colnames(mat)[df$col])

    print(df)
    if(!is.null(fn)) write.table(df, file = fn, quote = F, row.names = F, sep = "\t")
    return(df)
  }
