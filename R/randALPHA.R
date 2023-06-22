# ------------------------------------------------------------------------------------------
# Description: This function was lifted from the ebsRtools package version 0.2.0 which 
#              will be used to generate randomization for Alpha Lattice Design.
# Source: ebsRtools 0.2.0
# ------------------------------------------------------------------------------------------
#' @name randALPHA
#' @aliases randALPHA
#' @title randALPHA
#' 
#' @description Function that generates randomization for for Alpha Lattice Design.
#' @param trt	Number of treatments
#' @param k	Block size, that is: number of plots per block
#' @param r	Number of replications (blocks)
#' @param nTrial Number of occurrences
#' @param tag Format of plot number tag
#' @param continue logical
#' @source ebsRtools 0.2.0
#' 
# ------------------------------------------------------------------------------------------

randALPHA <- function (trt, k, r, tag = 2, nTrial = 1, continue = F) {
  number <- 10
  if (tag > 0) {
    number <- 10^tag
  }
  name.trt <- c(paste(deparse(substitute(trt))))
  v <- length(trt)
  if (!any(v%%2:(v - 1) == 0)) {
    stop("nTreatment is a prime number. The implemented Alpha-lattice designs (with equal block size) can not handle this")
  }
  else {
    s <- v/k
    if (v%%k != 0) 
      stop("\nThe block size is not appropriate. nTreatments must be multiple of k (block size).")
    else {
      for (w in c(1:nTrial)) {
        if (w == 1) {
          trials <- list()
        }
        E <- ((v - 1) * (r - 1))/((v - 1) * (r - 1) + 
                                    ((r) * (s - 1)))
        res.mod.s <- rep(c(0:s)[1:s], times = s)
        serie <- NA
        ex <- sample(res.mod.s, k * r, replace = T)
        if (r == 2 & k <= s) {
          serie <- "I"
          ex <- c(rep(0, k), seq(0, (k - 1)))
        }
        if (r == 3) {
          if (s%%2 != 0 & k <= s) {
            serie <- "II"
            ex <- c(rep(0, k), seq(0, (k - 1)), c(0, 
                                                  seq(max(res.mod.s), 0)[c(1:k - 1)]))
          }
          if (s%%2 == 0 & k <= s - 1) {
            serie <- "III"
            ex <- c(rep(0, k), seq(0, (k - 1)), as.vector(rbind(seq(0, 
                                                                    k)[seq(k/2)], seq(s/2, k * 2)[seq(k/2)])))
          }
        }
        if (r == 4 & s%%2 != 0 & s%%3 != 0 & k <= s) {
          serie <- "IV"
          ex <- c(rep(0, k), seq(0, (k - 1)), c(0, seq(s - 
                                                         1, 1))[c(1:(k))], as.vector(rbind(seq(0, 
                                                                                               (k + 1)), seq((s + 1)/2, k * 2)))[1:k])
        }
        alpha <- matrix(ex, k, r)
        alphax <- matrix(NA, k, s)
        for (o in c(1:r)) {
          for (l in c(1:s)) {
            alphax[, l] <- alpha[, o] + res.mod.s[l]
          }
          if (o == 1) {
            tmp <- alphax
          }
          else {
            tmp <- cbind(tmp, alphax)
          }
        }
        alphax <- tmp%%s
        plan <- alphax
        for (o in c(0:k - 1)) {
          plan[o + 1, ] <- alphax[o + 1, ] + (o * s)
        }
        if (v <= 1000) {
          M <- matrix(NA, v, v)
          rownames(M) <- c(1:v) - 1
          colnames(M) <- c(1:v) - 1
          for (i in as.numeric(rownames(M))) {
            for (j in as.numeric(colnames(M))) {
              if (j > i) {
                pi <- as.integer(i/s) + 1
                pj <- as.integer(j/s) + 1
                setconcur <- alpha[pj, ] - alpha[pi, 
                ]
                M[i + 1, j + 1] <- sum((j - i)%%s == 
                                         setconcur)
              }
            }
          }
          gs <- sort(unique(c(M))[!is.na(unique(c(M)))])
          alphaType <- paste("Alpha-Latice (", 
                             paste(gs, collapse = ","), ")", 
                             sep = "")
          report_gs <- matrix(NA, 1, length(gs))
          for (i in (1:length(gs))) {
            report_gs[i] <- sum(M == gs[i], na.rm = T)
          }
          Concurrences <- paste(gs, "(", report_gs, 
                                ")", sep = "")
        }
        else {
          alphaType <- "Alpha-Lattice"
          Concurrences <- "nTreatment too large, value not calculated"
        }
        book <- data.frame(origem = c(1:(r * v)), trt = as.numeric(plan) + 
                             1, block = rep(c(1:(r * s)), each = k), replicate = rep(c(1:r), 
                                                                                     each = v), randK = c(sample(1:(r * v), r * 
                                                                                                                   v)), randS = c(rep(sample(1:(s * r)), each = k)))
        book <- book[with(book, order(replicate, randS, 
                                      randK)), c("trt", "block", "replicate")]
        book$block <- rep(c(1:(r * s)), each = k)
        book$plot_number <- book$replicate * number + 
          (1:v)
        book$entry_id <- trt[book$trt]
        book <- book[, c("plot_number", "replicate", 
                         "block", "entry_id")]
        rownames(book) <- c(1:nrow(book))
        if (continue) {
          book$plot_number <- rownames(book)
        }
        parameters <- list(design = alphaType, Efficiency = E, 
                           trt = trt, k = k, r = r, Concurrences = Concurrences)
        design <- list(parameters = parameters, book = book)
        trials[[w]] <- design
      }
    }
  }
  return(trials)
}