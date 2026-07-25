## Synthetic .cool/.mcool fixtures.
##
## marinerData's .cool/.mcool example files are published on ExperimentHub
## (EH9511-EH9514, hosted at doi:10.5281/zenodo.10906240), but the released
## marinerData package builds its accessors from inst/extdata/metadata.csv,
## and those rows were never pushed to Bioconductor -- so the accessor
## functions do not exist in any installable version.
##
## Small local files are the better choice for tests regardless. They keep
## `R CMD check` off the network (the real .mcool is 55 MB) and, more usefully,
## let these tests assert on counts we chose ourselves rather than only
## checking that cooler output agrees with .hic output. Two real bugs in
## `coolStraw()` survived precisely because the original tests only did the
## latter.
##
## Layout follows the cooler spec: chroms/, bins/, pixels/, indexes/, with
## /resolutions/<binSize>/ nesting for .mcool.

.coolChroms <- function() {
    data.frame(name = c("1", "2"), length = c(1000L, 600L))
}

## Bin table for a given resolution. Bins are numbered globally across
## chromosomes, but `start`/`end` are chromosome-relative.
.coolBins <- function(binSize) {
    chroms <- .coolChroms()
    do.call(rbind, lapply(seq_len(nrow(chroms)), function(i) {
        n <- chroms$length[i] %/% binSize
        data.frame(
            chrom = rep(i - 1L, n),
            start = (seq_len(n) - 1L) * binSize,
            end   = seq_len(n) * binSize
        )
    }))
}

## 0-based bin offset of each chromosome, plus a terminator.
.coolChromOffset <- function(binSize) {
    as.integer(c(0L, cumsum(.coolChroms()$length %/% binSize)))
}

## Every upper-triangular bin pair, with a count encoding its own coordinates
## so a misplaced pixel is identifiable from its value alone.
.coolPixels <- function(nBins) {
    px <- do.call(rbind, lapply(0:(nBins - 1L), function(i) {
        j <- i:(nBins - 1L)
        data.frame(bin1_id = i, bin2_id = j, count = i * 1000L + j)
    }))
    px[order(px$bin1_id, px$bin2_id), ]
}

## 0-based pixel offset of each bin, plus a terminator.
.coolBin1Offset <- function(px, nBins) {
    as.integer(c(0L, cumsum(vapply(
        0:(nBins - 1L), function(i) sum(px$bin1_id == i), integer(1)
    ))))
}

## Non-uniform balancing weights, so BALANCE is distinguishable from NONE.
.coolWeights <- function(nBins) 1 / (seq_len(nBins) + 1)

.writeCoolGroups <- function(file, prefix, binSize) {
    g <- function(n) if (nzchar(prefix)) paste0(prefix, "/", n) else n

    chroms <- .coolChroms()
    bins <- .coolBins(binSize)
    nBins <- nrow(bins)
    px <- .coolPixels(nBins)

    rhdf5::h5createGroup(file, g("chroms"))
    rhdf5::h5write(chroms$name, file, g("chroms/name"))
    rhdf5::h5write(chroms$length, file, g("chroms/length"))

    rhdf5::h5createGroup(file, g("bins"))
    rhdf5::h5write(bins$chrom, file, g("bins/chrom"))
    rhdf5::h5write(bins$start, file, g("bins/start"))
    rhdf5::h5write(bins$end, file, g("bins/end"))
    rhdf5::h5write(.coolWeights(nBins), file, g("bins/weight"))

    rhdf5::h5createGroup(file, g("pixels"))
    rhdf5::h5write(px$bin1_id, file, g("pixels/bin1_id"))
    rhdf5::h5write(px$bin2_id, file, g("pixels/bin2_id"))
    rhdf5::h5write(px$count, file, g("pixels/count"))

    rhdf5::h5createGroup(file, g("indexes"))
    rhdf5::h5write(.coolChromOffset(binSize), file, g("indexes/chrom_offset"))
    rhdf5::h5write(.coolBin1Offset(px, nBins), file, g("indexes/bin1_offset"))
}

makeCoolFile <- function(binSize = 100L) {
    f <- tempfile(fileext = ".cool")
    rhdf5::h5createFile(f)
    .writeCoolGroups(f, "", binSize)
    rhdf5::h5closeAll()
    f
}

makeMcoolFile <- function(binSizes = c(100L, 200L)) {
    f <- tempfile(fileext = ".mcool")
    rhdf5::h5createFile(f)
    rhdf5::h5createGroup(f, "resolutions")
    for (bs in binSizes) {
        rhdf5::h5createGroup(f, paste0("resolutions/", bs))
        .writeCoolGroups(f, paste0("resolutions/", bs), bs)
    }
    rhdf5::h5closeAll()
    f
}

## Expected coolStraw() output for a query, derived independently of the
## implementation: pick the bins each range covers, take the upper-triangular
## pixels among them, and look up their chromosome-relative starts.
expectedCoolStraw <- function(binSize, chrom1, s1, e1, chrom2, s2, e2) {
    chroms <- .coolChroms()
    bins <- .coolBins(binSize)
    offs <- .coolChromOffset(binSize)
    px <- .coolPixels(nrow(bins))

    binsFor <- function(chrom, s, e) {
        i <- which(chroms$name == chrom)
        first <- offs[i]
        last <- offs[i + 1L] - 1L
        e <- min(e, chroms$length[i])
        seq(first + s %/% binSize, min(first + e %/% binSize, last))
    }

    keep <- px[px$bin1_id %in% binsFor(chrom1, s1, e1) &
                   px$bin2_id %in% binsFor(chrom2, s2, e2), ]
    out <- data.frame(
        x = bins$start[keep$bin1_id + 1L],
        y = bins$start[keep$bin2_id + 1L],
        counts = as.numeric(keep$count)
    )
    out[order(out$x, out$y), ]
}
