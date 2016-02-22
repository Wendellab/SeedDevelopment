## functions previously in WGCNA, not anymore since update to 1.4

plotMat<-
function (x, nrgcols = 50, rlabels = FALSE, clabels = FALSE, 
    rcols = 1, ccols = 1, title = "", ...) 
{
    n <- nrow(x)
    p <- ncol(x)
    image(1:p, 1:n, t(x[n:1, ]), col = rgcolors.func(nrgcols), 
        axes = FALSE, xlab = "", ylab = "", ...)
    if (length(ccols) == 1) {
        axis(3, at = 1:p, labels = clabels, las = 2, cex.axis = 0.6, 
            col.axis = ccols)
    }
    if (length(ccols) == p) {
        cols <- unique(ccols)
        for (i in 1:length(cols)) {
            which <- (1:p)[ccols == cols[i]]
            axis(3, at = which, labels = clabels[which], las = 2, 
                cex.axis = 0.6, col.axis = cols[i])
        }
    }
    if (length(rcols) == 1) {
        axis(2, at = n:1, labels = rlabels, las = 2, cex.axis = 0.6, 
            col.axis = rcols)
    }
    if (length(rcols) == n) {
        cols <- unique(rcols)
        for (i in 1:length(cols)) {
            which <- (1:n)[rcols == cols[i]]
            axis(2, at = (n:1)[which], labels = rlabels[which], 
                las = 2, cex.axis = 0.6, col.axis = cols[i])
        }
    }
    mtext(title, side = 3, line = 3)
    box()
}


rgcolors.func<-
function (n = 50) 
{
    k <- round(n/2)
    r <- c(rep(0, k), seq(0, 1, length = k))
    g <- c(rev(seq(0, 1, length = k)), rep(0, k))
    res <- rgb(r, g, rep(0, 2 * k))
    res
}

plotEigengeneNetworksCex<-function (multiME, setLabels, letterSubPlots = FALSE, Letters = NULL,
excludeGrey = TRUE, greyLabel = "grey", plotDendrograms = TRUE,
plotHeatmaps = TRUE, setMargins = TRUE, marDendro = NULL,
marHeatmap = NULL, colorLabels = TRUE, signed = TRUE, heatmapColors = NULL,
plotAdjacency = TRUE, printAdjacency = FALSE, cex.adjacency = 0.9,
coloredBarplot = TRUE, barplotMeans = TRUE, barplotErrors = FALSE,
plotPreservation = "standard", zlimPreservation = c(0, 1),
printPreservation = FALSE, cex.preservation = 0.9, cex=1,...)
{
    size = checkSets(multiME, checkStructure = TRUE)
    if (!size$structureOK) {
        multiME = fixDataStructure(multiME)
    }
    if (is.null(Letters))
    Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if (is.null(heatmapColors))
    if (signed) {
        heatmapColors = blueWhiteRed(50)
    }
    else {
        heatmapColors = heat.colors(30)
    }
    nSets = length(multiME)
    #cex = par("cex")
    mar = par("mar")
    nPlotCols = nSets
    nPlotRows = as.numeric(plotDendrograms) + nSets * as.numeric(plotHeatmaps)
    if (nPlotRows == 0)
    stop("Nothing to plot: neither dendrograms not heatmaps requested.")
    par(mfrow = c(nPlotRows, nPlotCols))
    par(cex = cex)
    if (excludeGrey)
    for (set in 1:nSets) multiME[[set]]$data = multiME[[set]]$data[,
    substring(names(multiME[[set]]$data), 3) != greyLabel]
    plotPresTypes = c("standard", "hyperbolic", "both")
    ipp = pmatch(plotPreservation, plotPresTypes)
    if (is.na(ipp))
    stop(paste("Invalid 'plotPreservation'. Available choices are",
    paste(plotPresTypes, sep = ", ")))
    letter.ind = 1
    if (plotDendrograms)
    for (set in 1:nSets) {
        par(mar = marDendro)
        labels = names(multiME[[set]]$data)
        uselabels = labels[substring(labels, 3) != greyLabel]
        corME = cor(multiME[[set]]$data[substring(labels,
        3) != greyLabel, substring(labels, 3) != greyLabel],
        use = "p")
        disME = as.dist(1 - corME)
        clust = fastcluster::hclust(disME, method = "average")
        if (letterSubPlots) {
            main = paste(substring(Letters, letter.ind, letter.ind),
            ". ", setLabels[set], sep = "")
        }
        else {
            main = setLabels[set]
        }
        plotLabels = uselabels
        plot(clust, main = main, sub = "", xlab = "", labels = plotLabels,
        ylab = "", ylim = c(0, 1))
        letter.ind = letter.ind + 1
    }
    if (plotHeatmaps)
    for (i.row in (1:nSets)) for (i.col in (1:nSets)) {
        letter.ind = i.row * nSets + i.col
        if (letterSubPlots) {
            letter = paste(substring(Letters, first = letter.ind,
            last = letter.ind), ".  ", sep = "")
        }
        else {
            letter = NULL
        }
        par(cex = cex)
        if (setMargins) {
            if (is.null(marHeatmap)) {
                if (colorLabels) {
                    par(mar = c(1, 2, 3, 4) + 0.2)
                }
                else {
                    par(mar = c(6, 7, 3, 5) + 0.2)
                }
            }
            else {
                par(mar = marHeatmap)
            }
        }
        nModules = dim(multiME[[i.col]]$data)[2]
        textMat = NULL
        if (i.row == i.col) {
            corME = cor(multiME[[i.col]]$data, use = "p")
            pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data))
            if (printAdjacency) {
                textMat = paste(signif(corME, 2), "\n", signif(pME,
                1))
                dim(textMat) = dim(corME)
            }
            if (signed) {
                if (plotAdjacency) {
                    if (printAdjacency) {
                        textMat = paste(signif((1 + corME)/2, 2),
                        "\n", signif(pME, 1))
                        dim(textMat) = dim(corME)
                    }
                    labeledHeatmap((1 + corME)/2, names(multiME[[i.col]]$data),
                    names(multiME[[i.col]]$data), main = paste(letter,
                    setLabels[[i.col]]), invertColors = FALSE,
                    zlim = c(0, 1), colorLabels = colorLabels,
                    colors = heatmapColors, setStdMargins = FALSE,
                    textMatrix = textMat, cex.text = cex.adjacency,
                    ...)
                }
                else {
                    labeledHeatmap(corME, names(multiME[[i.col]]$data),
                    names(multiME[[i.col]]$data), main = paste(letter,
                    setLabels[[i.col]]), invertColors = FALSE,
                    zlim = c(-1, 1), colorLabels = colorLabels,
                    colors = heatmapColors, setStdMargins = FALSE,
                    textMatrix = textMat, cex.text = cex.adjacency,
                    ...)
                }
            }
            else {
                labeledHeatmap(abs(corME), names(multiME[[i.col]]$data),
                names(multiME[[i.col]]$data), main = paste(letter,
                setLabels[[i.col]]), invertColors = FALSE,
                zlim = c(0, 1), colorLabels = colorLabels,
                colors = heatmapColors, setStdMargins = FALSE,
                textMatrix = textMat, cex.text = cex.adjacency,
                ...)
            }
        }
        else {
            corME1 = cor(multiME[[i.col]]$data, use = "p")
            corME2 = cor(multiME[[i.row]]$data, use = "p")
            cor.dif = (corME1 - corME2)/2
            d = tanh((corME1 - corME2)/(abs(corME1) + abs(corME2))^2)
            if (ipp == 1 | ipp == 3) {
                dispd = cor.dif
                main = paste(letter, "Preservation")
                if (ipp == 3) {
                    dispd[upper.tri(d)] = d[upper.tri(d)]
                    main = paste(letter, "Hyperbolic preservation (UT)\nStandard preservation (LT)")
                }
            }
            else {
                dispd = d
                main = paste(letter, "Hyperbolic preservation")
            }
            if (i.row > i.col) {
                if (signed) {
                    half = as.integer(length(heatmapColors)/2)
                    range = c(half:length(heatmapColors))
                    halfColors = heatmapColors[range]
                }
                else {
                    halfColors = heatmapColors
                }
                if (printPreservation) {
                    printMtx = matrix(paste(".", as.integer((1 -
                    abs(dispd)) * 100), sep = ""), nrow = nrow(dispd),
                    ncol = ncol(dispd))
                    printMtx[printMtx == ".100"] = "1"
                }
                else {
                    printMtx = NULL
                }
                if ((sum((1 - abs(dispd)) < zlimPreservation[1]) ||
                ((1 - abs(dispd)) > zlimPreservation[2])) >
                0)
                warning("plotEigengeneNetworks: Correlation preservation data out of zlim range.")
                labeledHeatmap(1 - abs(dispd), names(multiME[[i.col]]$data),
                names(multiME[[i.col]]$data), main = main,
                invertColors = FALSE, colorLabels = colorLabels,
                zlim = zlimPreservation, colors = halfColors,
                setStdMargins = FALSE, textMatrix = printMtx,
                cex.text = cex.preservation, ...)
            }
            else {
                if (ipp == 2) {
                    dp = 1 - abs(d)
                    method = "Hyperbolic:"
                }
                else {
                    dp = 1 - abs(cor.dif)
                    method = "Preservation:"
                }
                diag(dp) = 0
                if (barplotMeans) {
                    sum_dp = mean(dp[upper.tri(dp)])
                    means = apply(dp, 2, sum)/(ncol(dp) - 1)
                    if (barplotErrors) {
                        errors = sqrt((apply(dp^2, 2, sum)/(ncol(dp) -
                        1) - means^2)/(ncol(dp) - 2))
                    }
                    else {
                        errors = NULL
                    }
                    labeledBarplot(means, names(multiME[[i.col]]$data),
                    main = paste(letter, "D=", signif(sum_dp,
                    2)), ylim = c(0, 1), colorLabels = colorLabels,
                    colored = coloredBarplot, setStdMargins = FALSE,
                    stdErrors = errors, ...)
                }
                else {
                    sum_dp = sum(dp[upper.tri(dp)])
                    labeledBarplot(dp, names(multiME[[i.col]]$data),
                    main = paste(letter, method, "sum = ",
                    signif(sum_dp, 3)), ylim = c(0, dim(dp)[[1]]),
                    colorLabels = colorLabels, colored = coloredBarplot,
                    setStdMargins = FALSE, ...)
                }
            }
        }
    }
}
