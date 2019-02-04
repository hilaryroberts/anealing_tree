
entropy <- function(vec) {
    ratios <- table(vec)/length(vec)
    out <- log(ratios, base = length(ratios))*ratios
    out[is.nan(out)] <- 0
    -sum(out)
}

ig <- function(bins, response){
    weights <- table(bins)/length(bins)
    children <- tapply(response, INDEX = bins, FUN = entropy)
    entropy(response) - sum(weights * children)
}


find_split <- function(var, response) {
    var <- sort(var)
    response <- response[order(var)]
    out <- c()
    for (i in seq_along(var)) {
        bins <- var < var[i]
        out <- append(out, ig(bins, response))
    }
    split <- list()
    split[['point']] <- var[out == max(out)][1]
    split[['ig']] <- max(out)
    split
}

factor_split <- function(var, response){
    ranking <- sapply(levels(var), function(l) {
        bins <- var == l
        ig(bins, response)
    })
    split <- list()
    split[['point']] <- names(which(ranking == max(ranking)))[1]
    split[['ig']] <- max(ranking)
    split
}


rank_vars <- function(predictors, response){
    splits <- c()
    for (v in colnames(predictors)) {
        if (is.numeric(predictors[[v]])) {
            split <- find_split(predictors[[v]], response)
            predictors[[v]] <- as.factor(ifelse(predictors[[v]] < split$point, paste0('<', split$point), paste0('>=', split$point)))
        } else if (is.factor(predictors[[v]])) {
            split <- factor_split(predictors[[v]], response)
            predictors[[v]] <- as.factor(ifelse(predictors[[v]] == split$point, split$point, paste('not', split$point)))
        } else stop(paste(v, 'has unrecognised variable type'))
        splits <- c(splits, list(split))
    }
    names(splits) <- colnames(predictors)
    ranking <- sapply(splits, FUN = function(x) x$ig)
    ordering <- order(ranking, decreasing = TRUE)
    splits[ordering]
}


#replace with method?

split_node <- function(data, respname, n = 1, igthresh = 0){
    ranking <- rank_vars(data[, -which(colnames(data) == respname)], data[[respname]])
    split <- ranking[[n]]['point']
    var <- data[[names(ranking)[n]]]
    if (ranking[[n]]['ig'] <= igthresh) {
        return('Close')
    }
    if (is.numeric(var)) {
        lchild <- data[var < as.numeric(split), ]
        rchild <- data[var >= as.numeric(split), ]
        out <- list(lchild, rchild)
        names(out) <- c(paste0(names(ranking[n]), '<', split), paste0(names(ranking[n]), '>=', split))
    } else {
        lchild <- data[var == split, ]
        rchild <- data[var != split, ]
        out <- list(lchild, rchild)
        names(out) <- c(paste0(names(ranking[n]), '=="', split, '"'), paste0(names(ranking[n]), '!="', split, '"'))
    }
    out
}


library(data.tree)

build_tree <- function(x, respname, maxdepth = Inf, minig = 0){
    outtree <- Node$new('rootnode')
    outtree$data <- x
    xnode <- outtree
    open <- TRUE
    while (any(open)) {
        newnodes <- split_node(xnode$data, respname = respname, igthresh = minig)
        if (newnodes == 'Close' | xnode$level > maxdepth) {
            xnode$closed <- TRUE
        } else {
            xnode$AddChild(names(newnodes)[1], data = newnodes[[1]])
            xnode$AddChild(names(newnodes)[2], data = newnodes[[2]])
            if (sum(table(xnode$children[1][[1]]$data[[respname]]) != 0) < 2) xnode$children[1][[1]]$closed <- TRUE
            if (sum(table(xnode$children[2][[1]]$data[[respname]]) != 0) < 2) xnode$children[2][[1]]$closed <- TRUE
        }
        open <- sapply(outtree$leaves, function(x) is.null(x$closed))
        xnode <- outtree$leaves[open][1][[1]]
    }
    for (l in outtree$leaves) {
        dis <- table(l$data[[respname]])
        l$decision <- names(dis)[dis == max(dis)]
    }
    outtree$Do(function(node) node$data <- NULL)
    outtree
}


aneal <- function(l, n) {
    throws <- runif(n)
    threshs <- l/(n:1)
    which(throws < threshs)[1]
}



#DROP LEVELS WITH 0 obs

anealing_tree <- function(x, respname, maxdepth = Inf, minig = 0){
    outtree <- Node$new('rootnode')
    outtree$data <- x
    xnode <- outtree
    open <- TRUE
    nvar <- ncol(x) - 1
    while (any(open)) {
        af <- aneal(xnode$level, nvar)
        newnodes <- split_node(xnode$data, respname = respname, n = af, igthresh = minig)
        if (newnodes == 'Close' | xnode$level > maxdepth) {
            xnode$closed <- TRUE
        } else {
            xnode$AddChild(names(newnodes)[1], data = newnodes[[1]])
            xnode$AddChild(names(newnodes)[2], data = newnodes[[2]])
            if (sum(table(xnode$children[1][[1]]$data[[respname]]) != 0) < 2) xnode$children[1][[1]]$closed <- TRUE
            if (sum(table(xnode$children[2][[1]]$data[[respname]]) != 0) < 2) xnode$children[2][[1]]$closed <- TRUE
        }
        open <- sapply(outtree$leaves, function(x) is.null(x$closed))
        xnode <- outtree$leaves[open][1][[1]]
    }
    for (l in outtree$leaves) {
        dis <- table(l$data[[respname]])
        l$decision <- names(dis)[dis == max(dis)]
    }
    outtree$Do(function(node) node$data <- NULL)
    outtree
}



tree_predict <- function(x, model){
    xnode <- model$root
    while (!xnode$isLeaf) {
        condition <- xnode$children[[1]]$name
        if (with(x, eval(parse(text = condition)))) {
            xnode <- xnode$children[[1]]
        } else xnode <- xnode$children[[2]]
    }
    xnode$decision
}

multi_predict <- function(model, x){
    rowlist <- split(x, 1:nrow(x))
    sapply(rowlist, FUN = tree_predict, model = model)
}

classic_forest <- function(data, respname, ntree, mtry, minig = 0, maxdepth = Inf) {
    forest <- list()
    predictors <- colnames(data)[colnames(data) != respname]
    for (i in 1:ntree) {
        samp <- predictors[sample(length(predictors), mtry)]
        sampdata <- data[, c(respname, samp)]
        tree <- build_tree(sampdata, respname, minig = minig, maxdepth = maxdepth)
        forest <- append(forest, tree)
    }
    forest
}


anealing_forest <- function(data, respname, ntree, minig = 0, maxdepth = Inf) {
    forest <- list()
    for (i in 1:ntree) {
        tree <- anealing_tree(data, respname, minig = minig, maxdepth = maxdepth)
        forest <- append(forest, tree)
    }
    forest
}

forest_predict <- function(x, model, result) {
    votelist <- lapply(model, multi_predict, x)
    votes <- do.call('rbind', votelist)
    #tally <- lapply(as.data.frame(votes), table)
    #majorities <- sapply(tally, function(x) names(x)[x == max(x)][1])
    scores <- sapply(as.data.frame(votes), function(x) mean(x == result))
    scores
}
