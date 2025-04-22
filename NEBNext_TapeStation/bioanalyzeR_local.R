# Local implementation of bioanalyzeR functions
# Auto-generated - Do not modify manually

# Required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
   install.packages("ggplot2")
 }
if (!requireNamespace("gridExtra", quietly = TRUE)) {
   install.packages("gridExtra")
 }
if (!requireNamespace("grid", quietly = TRUE)) {
   install.packages("grid")
 }
if (!requireNamespace("reshape2", quietly = TRUE)) {
   install.packages("reshape2")
 }
if (!requireNamespace("dplyr", quietly = TRUE)) {
   install.packages("dplyr")
 }
if (!requireNamespace("plyr", quietly = TRUE)) {
   install.packages("plyr")
 }
if (!requireNamespace("scales", quietly = TRUE)) {
   install.packages("scales")
 }
if (!requireNamespace("XML", quietly = TRUE)) {
   install.packages("XML")
 }
if (!requireNamespace("methods", quietly = TRUE)) {
   install.packages("methods")
 }

# Load required packages
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(dplyr)
library(plyr)
library(scales)
library(XML)
library(methods)

# bioanalyzeR internal functions

# Extracted functions

LOWER.MARKER.NAMES <-
  c("Lower Marker", "edited Lower Marker")
UPPER.MARKER.NAMES <-
  c("Upper Marker", "edited Upper Marker")
NROW.ALIGNED <-
  760

# Function: read.tapestation
read.tapestation <- function (xml.file, csv.file = NULL, method = "hyman", extrapolate = FALSE) 
{
    if (is.null(csv.file)) {
        candidate.files <- sapply(c(csv1 = "_Electropherogram.csv", 
            csv2 = "_Electropherogram.csv.gz"), function(suffix) sub("\\.xml(\\.gz)?$", 
            suffix, xml.file))
        found.files <- sapply(candidate.files, file.exists)
        stopifnot(`found conflicting CSV files` = sum(found.files) <= 
            1)
        stopifnot(`couldn't find CSV file` = sum(found.files) > 
            0)
        csv.file <- candidate.files[found.files]
    }
    parsed.data <- read.tapestation.xml(xml.file)
    stopifnot(`multiple batches provided` = length(unique(parsed.data$samples$batch)) == 
        1)
    batch <- parsed.data$samples$batch[1]
    result <- electrophoresis(data = read.tapestation.csv(csv.file), 
        assay.info = setNames(list(parsed.data$assay.info), batch), 
        samples = parsed.data$samples, peaks = parsed.data$peaks, 
        regions = parsed.data$regions)
    if (sum(result$samples$is.ladder) > 1 && length(unique(result$samples$ladder.well)) == 
        1) {
        result <- subset(result, !(is.ladder & well.number != 
            ladder.well))
    }
    is.lower.marker <- result$peaks$peak.observations %in% LOWER.MARKER.NAMES
    marker.distances <- data.frame(lower = sapply(seq(nrow(result$samples)), 
        function(i) {
            distance <- result$peaks$distance[is.lower.marker & 
                result$peaks$sample.index == i]
            if (length(distance) == 0) {
                return(NA)
            }
            else if (length(distance) == 1) {
                return(distance)
            }
            else {
                stop(paste("multiple lower marker peaks for sample", 
                  i))
            }
        }))
    is.upper.marker <- result$peaks$peak.observations %in% UPPER.MARKER.NAMES
    if (sum(is.upper.marker) == 0) {
        marker.distances$upper <- 0
    }
    else {
        marker.distances$upper <- sapply(seq(nrow(result$samples)), 
            function(i) {
                distance <- result$peaks$distance[is.upper.marker & 
                  result$peaks$sample.index == i]
                if (length(distance) == 0) {
                  return(NA)
                }
                else if (length(distance) == 1) {
                  return(distance)
                }
                else {
                  stop(paste("multiple upper marker peaks for sample", 
                    i))
                }
            })
    }
    marker.distances$range <- marker.distances$lower - marker.distances$upper
    result$data$relative.distance <- (result$data$distance - 
        marker.distances$upper[result$data$sample.index])/marker.distances$range[result$data$sample.index]
    result$peaks$relative.distance <- (result$peaks$distance - 
        marker.distances$upper[result$peaks$sample.index])/marker.distances$range[result$peaks$sample.index]
    result$peaks$lower.relative.distance <- (result$peaks$lower.distance - 
        marker.distances$upper[result$peaks$sample.index])/marker.distances$range[result$peaks$sample.index]
    result$peaks$upper.relative.distance <- (result$peaks$upper.distance - 
        marker.distances$upper[result$peaks$sample.index])/marker.distances$range[result$peaks$sample.index]
    if (sum(!is.na(result$samples$ladder.well)) == 0) {
        warning("no ladder identified so calibrations are not performed")
        return(result)
    }
    result <- calculate.molarity(calculate.concentration(calculate.length(result, 
        method, extrapolate)))
    if (!is.null(result$regions)) {
        result$regions$lower.distance <- result$regions$lower.relative.distance * 
            marker.distances$range[result$regions$sample.index] + 
            marker.distances$upper[result$regions$sample.index]
        result$regions$upper.distance <- result$regions$upper.relative.distance * 
            marker.distances$range[result$regions$sample.index] + 
            marker.distances$upper[result$regions$sample.index]
    }
    result$samples$is.ladder <- NULL
    result
}

# Function: calculate.concentration
calculate.concentration <- function (electrophoresis, ladder.concentrations = NULL) 
{
    x.name <- get.x.name(electrophoresis)
    delta <- do.call(rbind, by(electrophoresis$data, electrophoresis$data$sample.index, 
        function(data.subset) data.frame(fluorescence = c(NA, 
            diff(data.subset$fluorescence)), x = c(NA, diff(data.subset[[x.name]]))), 
        simplify = F))
    if (x.name == "relative.distance") 
        delta$x <- -delta$x
    electrophoresis$data$area <- (2 * electrophoresis$data$fluorescence - 
        delta$fluorescence) * delta$x
    if (x.name == "aligned.time") 
        electrophoresis$data$area <- electrophoresis$data$area/electrophoresis$data[[x.name]]
    has.upper.marker <- any(electrophoresis$peaks$peak.observations %in% 
        UPPER.MARKER.NAMES)
    which.markers <- lapply(seq(nrow(electrophoresis$samples)), 
        function(sample.index) which(electrophoresis$peaks$sample.index == 
            sample.index & ((electrophoresis$peaks$peak.observations %in% 
            LOWER.MARKER.NAMES & (if (is.null(ladder.concentrations)) 
            TRUE
        else electrophoresis$peaks$concentration == ladder.concentrations[1])) | 
            (electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES & 
                (if (is.null(ladder.concentrations)) 
                  TRUE
                else electrophoresis$peaks$concentration == ladder.concentrations[length(ladder.concentrations)])))))
    for (i in seq(nrow(electrophoresis$samples))) if (length(which.markers[[i]]) == 
        0) 
        warning(paste("no markers detected in sample", i))
    mass.coefficients <- rep(NA, nrow(electrophoresis$samples))
    for (batch in unique(electrophoresis$samples$batch)) {
        in.this.batch <- electrophoresis$samples$batch == batch
        for (ladder.well in unique(electrophoresis$samples$ladder.well[which(in.this.batch)])) {
            ladder.index <- which(in.this.batch & electrophoresis$samples$well.number == 
                ladder.well)
            stopifnot(`conflicting ladders` = length(ladder.index) == 
                1)
            which.samples <- which(in.this.batch & electrophoresis$samples$ladder.well == 
                ladder.well)
            this.ladder.concentrations <- if (!is.null(ladder.concentrations)) 
                ladder.concentrations
            else subset(electrophoresis$peaks, sample.index == 
                ladder.index)$concentration
            non.marker.concentrations <- this.ladder.concentrations[-c(1, 
                if (has.upper.marker) length(this.ladder.concentrations) else NULL)]
            ladder.peaks <- which(electrophoresis$peaks$sample.index == 
                ladder.index & (if (is.null(ladder.concentrations)) 
                TRUE
            else electrophoresis$peaks$concentration %in% ladder.concentrations))
            stopifnot(`no ladder peaks recognized` = length(ladder.peaks) > 
                0)
            non.marker.peaks <- setdiff(ladder.peaks, which.markers[[ladder.index]])
            ladder.areas <- integrate.peak(electrophoresis, ladder.peaks, 
                "area")
            non.marker.areas <- ladder.areas[!ladder.peaks %in% 
                which.markers[[ladder.index]]]
            ladder.marker.areas <- ladder.areas[ladder.peaks %in% 
                which.markers[[ladder.index]]]
            ladder.mass.coefficient <- mean(electrophoresis$peaks$concentration[non.marker.peaks]/non.marker.areas, 
                na.rm = T)
            for (sample.index in which.samples) {
                sample.marker.areas <- integrate.peak(electrophoresis, 
                  which.markers[[sample.index]], "area")
                mass.coefficients[sample.index] <- if (length(sample.marker.areas) == 
                  0) 
                  NA
                else ladder.mass.coefficient * mean(ladder.marker.areas/sample.marker.areas, 
                  na.rm = T)
            }
        }
    }
    electrophoresis$data$concentration <- mass.coefficients[electrophoresis$data$sample.index] * 
        electrophoresis$data$area
    electrophoresis$data$area <- NULL
    electrophoresis
}

# Function: electrophoresis
electrophoresis <- function (data = NULL, assay.info = NULL, samples = NULL, peaks = NULL, 
    regions = NULL, calibration = NULL) 
structure(list(data = data, assay.info = assay.info, samples = samples, 
    peaks = peaks, regions = regions, calibration = calibration), 
    class = "electrophoresis")

# Function: get.x.name
get.x.name <- function (electrophoresis, raw = FALSE, allow.multiple = FALSE) 
{
    possible.x.names <- if (raw) 
        c("time", "distance")
    else c("aligned.time", "relative.distance")
    result <- intersect(possible.x.names, colnames(electrophoresis$data))
    stopifnot(`no x-variable found` = length(result) > 0, `multiple x-variables` = allow.multiple || 
        length(result) == 1)
    result
}

# Function: integrate.peak
integrate.peak <- function (electrophoresis, index = seq(nrow(electrophoresis$peaks)), 
    sum.variable = "molarity") 
sapply(index, function(i) sum(electrophoresis$data[[sum.variable]][in.peak(electrophoresis, 
    i)]))

# Function: in.peak
in.peak <- function (electrophoresis, which.peak) 
{
    x.names <- get.x.name(electrophoresis, allow.multiple = T)
    peak.x.name <- x.names[which(!is.na(electrophoresis$peaks[which.peak, 
        x.names]))]
    if (length(peak.x.name) == 0) {
        rep(NA, nrow(electrophoresis$data))
    }
    else {
        stopifnot(`multiple peak x-values` = length(peak.x.name) == 
            1)
        electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[which.peak] & 
            electrophoresis$data[[peak.x.name]] >= electrophoresis$peaks[[paste0("lower.", 
                peak.x.name)]][which.peak] & electrophoresis$data[[peak.x.name]] <= 
            electrophoresis$peaks[[paste0("upper.", peak.x.name)]][which.peak]
    }
}

# Function: subset.electrophoresis
subset.electrophoresis <- function (electrophoresis, ...) 
{
    nrow.initial <- nrow(electrophoresis$samples)
    electrophoresis$samples <- subset(electrophoresis$samples, 
        ...)
    if (nrow(electrophoresis$samples) == nrow.initial) {
        return(electrophoresis)
    }
    else if (nrow(electrophoresis$samples) == 0) {
        for (member in names(electrophoresis)) electrophoresis[[member]] <- NULL
        return(electrophoresis)
    }
    remaining.samples <- as.integer(rownames(electrophoresis$samples))
    electrophoresis$data <- subset(electrophoresis$data, sample.index %in% 
        remaining.samples)
    electrophoresis$assay.info <- electrophoresis$assay.info[names(electrophoresis$assay.info) %in% 
        as.character(electrophoresis$samples$batch)]
    if (!is.null(electrophoresis$peaks)) {
        electrophoresis$peaks <- subset(electrophoresis$peaks, 
            sample.index %in% remaining.samples)
        if (nrow(electrophoresis$peaks) == 0) 
            electrophoresis$peaks <- NULL
    }
    if (!is.null(electrophoresis$regions)) {
        electrophoresis$regions <- subset(electrophoresis$regions, 
            sample.index %in% remaining.samples)
        if (nrow(electrophoresis$regions) == 0) 
            electrophoresis$regions <- NULL
    }
    new.indices <- setNames(seq(nrow(electrophoresis$samples)), 
        rownames(electrophoresis$samples))
    electrophoresis$data$sample.index <- new.indices[as.character(electrophoresis$data$sample.index)]
    if (!is.null(electrophoresis$peaks)) 
        electrophoresis$peaks$sample.index <- new.indices[as.character(electrophoresis$peaks$sample.index)]
    if (!is.null(electrophoresis$regions)) 
        electrophoresis$regions$sample.index <- new.indices[as.character(electrophoresis$regions$sample.index)]
    rownames(electrophoresis$data) <- seq(nrow(electrophoresis$data))
    rownames(electrophoresis$samples) <- new.indices
    if (!is.null(electrophoresis$peaks)) 
        rownames(electrophoresis$peaks) <- seq(nrow(electrophoresis$peaks))
    if (!is.null(electrophoresis$regions)) 
        rownames(electrophoresis$regions) <- seq(nrow(electrophoresis$regions))
    electrophoresis
}

# Function: calculate.length
calculate.length <- function (electrophoresis, method = union(c("hyman", "interpolation", 
    "loglinear"), eval(formals(splinefun)$method)), extrapolate = FALSE) 
{
    method <- match.arg(method)
    x.name <- get.x.name(electrophoresis)
    lower.name <- paste0("lower.", x.name)
    upper.name <- paste0("upper.", x.name)
    electrophoresis$data$length <- NA
    electrophoresis$peaks$lower.length <- NA
    electrophoresis$peaks$upper.length <- NA
    if (!is.null(electrophoresis$regions)) {
        electrophoresis$regions[[lower.name]] <- NA
        electrophoresis$regions[[upper.name]] <- NA
    }
    electrophoresis$calibration <- list()
    for (batch in unique(electrophoresis$samples$batch)) {
        electrophoresis$calibration[[batch]] <- list()
        in.this.batch <- electrophoresis$samples$batch == batch
        for (ladder.well in unique(electrophoresis$samples$ladder.well[which(in.this.batch)])) {
            which.ladder.index <- which(in.this.batch & electrophoresis$samples$well.number == 
                ladder.well)
            peaks.ladder <- subset(electrophoresis$peaks, sample.index == 
                which.ladder.index)
            which.samples <- which(in.this.batch & electrophoresis$samples$ladder.well == 
                ladder.well)
            which.rows <- which(electrophoresis$data$sample.index %in% 
                which.samples)
            which.peaks <- which(electrophoresis$peaks$sample.index %in% 
                which.samples)
            which.regions <- which(electrophoresis$regions$sample.index %in% 
                which.samples)
            peaks.ladder$x <- peaks.ladder[[x.name]]
            if (method == "interpolation") {
                warning("linear interpolation gives ugly results for molarity estimation")
                standard.curve.function <- approxfun(peaks.ladder$x, 
                  peaks.ladder$length)
                standard.curve.inverse <- approxfun(peaks.ladder$length, 
                  peaks.ladder$x)
            }
            else if (method == "loglinear") {
                if (x.name == "relative.distance") {
                  mobility.model <- lm(x ~ log(length), peaks.ladder)
                  coefs <- coefficients(mobility.model)
                  standard.curve.function <- function(relative.distance) exp((relative.distance - 
                    coefs[1])/coefs[2])
                  standard.curve.inverse <- function(length) coefs[1] + 
                    coefs[2] * log(length)
                }
                else if (x.name == "aligned.time") {
                  mobility.model <- lm(1/aligned.time ~ log(length), 
                    data = peaks.ladder)
                  coefs <- coefficients(mobility.model)
                  standard.curve.function <- function(aligned.time) exp((1/aligned.time - 
                    coefs[1])/coefs[2])
                  standard.curve.inverse <- function(length) 1/(coefs[1] + 
                    log(length) * coefs[2])
                }
            }
            else {
                standard.curve.function <- splinefun(peaks.ladder$x, 
                  peaks.ladder$length, method = method)
                standard.curve.inverse <- splinefun(peaks.ladder$length, 
                  peaks.ladder$x, method = method)
            }
            electrophoresis$calibration[[batch]][[ladder.well]] <- list(mobility.function = standard.curve.function, 
                mobility.inverse = standard.curve.inverse, ladder.peaks = setNames(data.frame(peaks.ladder$length, 
                  peaks.ladder$x), c("length", x.name)))
            electrophoresis$data$length[which.rows] <- standard.curve.function(electrophoresis$data[[x.name]][which.rows])
            if (!extrapolate) {
                electrophoresis$data$length[!(in.custom.region(electrophoresis$data, 
                  min(peaks.ladder$length), max(peaks.ladder$length)) & 
                  in.custom.region(electrophoresis$data, min(peaks.ladder$x), 
                    max(peaks.ladder$x), bound.variable = x.name))] <- NA
            }
            else {
                if (method != "loglinear") 
                  warning("mobility model was fit by interpolation so there is no expectation of accuracy outside the ladder range")
                for (i in unique(electrophoresis$data$sample.index)) {
                  which.negative <- which(electrophoresis$data$sample.index == 
                    i & electrophoresis$data$length < 0)
                  if (any(which.negative, na.rm = T)) 
                    electrophoresis$data$length[seq(which(electrophoresis$data$sample.index == 
                      i)[1], tail(which.negative, 1))] <- NA
                }
            }
            if (x.name == "relative.distance") {
                lower.length.analog <- upper.name
                upper.length.analog <- lower.name
            }
            else if (x.name == "aligned.time") {
                lower.length.analog <- lower.name
                upper.length.analog <- upper.name
            }
            electrophoresis$peaks$lower.length[which.peaks] <- standard.curve.function(electrophoresis$peaks[[lower.length.analog]][which.peaks])
            electrophoresis$peaks$upper.length[which.peaks] <- standard.curve.function(electrophoresis$peaks[[upper.length.analog]][which.peaks])
            if (!is.null(electrophoresis$regions)) {
                electrophoresis$regions[[lower.length.analog]][which.regions] <- standard.curve.inverse(electrophoresis$regions$lower.length[which.regions])
                electrophoresis$regions[[upper.length.analog]][which.regions] <- standard.curve.inverse(electrophoresis$regions$upper.length[which.regions])
            }
        }
        electrophoresis$assay.info[[batch]]$method <- method
    }
    electrophoresis
}

# Function: in.custom.region
in.custom.region <- function (data, lower.bound = -Inf, upper.bound = Inf, bound.variable = "length") 
!is.na(data[[bound.variable]]) & data[[bound.variable]] >= lower.bound & 
    data[[bound.variable]] <= upper.bound

# Function: calculate.molarity
calculate.molarity <- function (electrophoresis) 
{
    electrophoresis$data$molarity <- NA
    for (batch in unique(electrophoresis$samples$batch)) {
        which.rows <- which(electrophoresis$data$sample.index %in% 
            which(electrophoresis$samples$batch == batch))
        electrophoresis$data$molarity[which.rows] <- electrophoresis$data$concentration[which.rows]/molecular.weight(electrophoresis$data$length[which.rows], 
            electrophoresis$assay.info[[batch]]$assay.type) * 
            1e+06
    }
    electrophoresis
}

# Function: molecular.weight
molecular.weight <- function (length, type) 
switch(type, DNA = length * 607.4 + 157.9, RNA = length * 320.5 + 
    159)

# Function: read.tapestation.csv
read.tapestation.csv <- function (csv.file) 
{
    raw.data <- read.csv(csv.file, encoding = "latin1")
    if (nrow(raw.data) == NROW.ALIGNED) 
        warning(paste("this looks like \"aligned\" data; did you export the unaligned data correctly?", 
            csv.file))
    data.frame(sample.index = rep(seq(ncol(raw.data)), each = nrow(raw.data)), 
        distance = rep(rev(seq(nrow(raw.data))), ncol(raw.data))/nrow(raw.data), 
        fluorescence = unlist(raw.data), row.names = seq(prod(dim(raw.data))))
}

# Function: read.tapestation.xml
read.tapestation.xml <- function (xml.file) 
{
    batch <- sub("\\.xml(\\.gz)?$", "", basename(xml.file))
    xml.root <- xmlRoot(xmlParse(xml.file))
    assay.info <- list(file.name = xmlValue(xml.root[["FileInformation"]][["FileName"]]), 
        creation.date = xmlValue(xml.root[["FileInformation"]][["RunEndDate"]]), 
        assay.name = xmlValue(xml.root[["FileInformation"]][["Assay"]]), 
        assay.type = NULL, length.unit = NULL, concentration.unit = NULL, 
        molarity.unit = NULL)
    if (grepl("RNA", assay.info$assay.name)) {
        assay.info$assay.type <- "RNA"
    }
    else if (grepl("D", assay.info$assay.name)) {
        assay.info$assay.type <- "DNA"
    }
    else {
        stop("unrecognized assay name")
    }
    assay.info$length.unit <- xmlValue(xml.root[["Assay"]][["Units"]][["MolecularWeightUnit"]])
    assay.info$concentration.unit <- xmlValue(xml.root[["Assay"]][["Units"]][["ConcentrationUnit"]])
    assay.info$molarity.unit <- switch(assay.info$concentration.unit, 
        `ng/µl` = "nM", `pg/µl` = "pM")
    result.list <- xmlApply(xml.root[["Samples"]], function(sample.xml) {
        well.number <- xmlValue(sample.xml[["WellNumber"]])
        sample.name <- trimws(xmlValue(sample.xml[["Comment"]]))
        sample.observations <- trimws(xmlValue(sample.xml[["Observations"]]))
        if (sample.observations == "Marker(s) not detected") 
            warning(paste(sample.observations, "for well", well.number, 
                sample.name))
        is.ladder <- grepl("Ladder(?! run as sample)(?! sizing changed)", 
            sample.observations, perl = T)
        if (sample.name == "") 
            sample.name <- well.number
        well.row <- ifelse(is.ladder && sample.name == "Electronic Ladder", 
            NA, substr(well.number, 1, 1))
        well.col <- ifelse(is.ladder && sample.name == "Electronic Ladder", 
            NA, substr(well.number, 2, 3))
        suppressWarnings(RINe <- as.numeric(xmlValue(sample.xml[["RNA"]][["RINe"]])))
        suppressWarnings(DIN <- as.numeric(xmlValue(sample.xml[["DIN"]])))
        reagent.id <- xmlValue(sample.xml[["ScreenTapeID"]])
        suppressWarnings(peaks <- if (length(xmlChildren(sample.xml[["Peaks"]])) == 
            0) 
            NULL
        else {
            peaks.raw <- xmlToDataFrame(sample.xml[["Peaks"]], 
                stringsAsFactors = F)
            data.frame(peak.observations = trimws(peaks.raw$Observations), 
                peak.comment = trimws(peaks.raw$Comment), length = as.integer(peaks.raw$Size), 
                distance = as.numeric(peaks.raw$RunDistance)/100, 
                lower.distance = as.numeric(peaks.raw$FromPercent)/100, 
                upper.distance = as.numeric(peaks.raw$ToPercent)/100, 
                concentration = as.numeric(peaks.raw$CalibratedQuantity), 
                molarity = as.numeric(peaks.raw$Molarity), stringsAsFactors = F)
        })
        suppressWarnings(regions <- if (length(xmlChildren(sample.xml[["Regions"]])) == 
            0) 
            NULL
        else {
            regions.raw <- xmlToDataFrame(sample.xml[["Regions"]], 
                stringsAsFactors = F)
            data.frame(region.comment = trimws(regions.raw$Comment), 
                lower.length = as.integer(regions.raw$From), 
                upper.length = as.integer(regions.raw$To), average.length = as.integer(regions.raw$AverageSize), 
                concentration = as.integer(regions.raw$Concentration), 
                molarity = as.numeric(regions.raw$Molarity), 
                proportion.of.total = as.numeric(regions.raw$PercentOfTotal)/100, 
                stringsAsFactors = F)
        })
        electrophoresis(samples = data.frame(batch, well.number, 
            well.row, well.col, sample.name, sample.observations, 
            reagent.id, RINe, DIN, is.ladder, stringsAsFactors = F), 
            peaks = peaks, regions = regions)
    })
    has.peaks <- !all(sapply(result.list, function(x) is.null(x$peaks)))
    has.regions <- !all(sapply(result.list, function(x) is.null(x$regions)))
    result <- electrophoresis(samples = do.call(rbind, c(lapply(result.list, 
        function(x) x$samples), make.row.names = F)), peaks = if (!has.peaks) 
        NULL
    else do.call(rbind, c(lapply(seq_along(result.list), function(i) if (is.null(result.list[[i]]$peaks)) NULL else cbind(sample.index = i, 
        result.list[[i]]$peaks)), make.row.names = F)), regions = if (!has.regions) 
        NULL
    else do.call(rbind, c(lapply(seq_along(result.list), function(i) if (is.null(result.list[[i]]$regions)) NULL else cbind(sample.index = i, 
        result.list[[i]]$regions)), make.row.names = F)), assay.info = assay.info)
    if (all(is.na(result$samples$RINe))) 
        result$samples$RINe <- NULL
    if (all(is.na(result$samples$DIN))) 
        result$samples$DIN <- NULL
    for (field in c("batch", "well.number", "sample.name", "reagent.id", 
        "sample.observations")) result$samples[, field] <- factor(result$samples[, 
        field], levels = unique(result$samples[, field]))
    result$samples$well.row <- factor(result$samples$well.row, 
        levels = LETTERS[1:8])
    result$samples$well.col <- factor(result$samples$well.col, 
        levels = 1:12)
    if (has.peaks) 
        for (field in c("peak.observations", "peak.comment")) result$peaks[, 
            field] <- factor(result$peaks[, field])
    if (has.regions) 
        for (field in c("region.comment")) result$regions[, field] <- factor(result$regions[, 
            field])
    ladder.wells <- result$samples$well.number[which(result$samples$is.ladder)]
    result$samples$ladder.well <- factor(NA, levels = levels(result$samples$well.number))
    if (length(ladder.wells) == 0) {
        return(result)
    }
    else if (length(ladder.wells) == 1) {
        result$samples$ladder.well <- ladder.wells
    }
    else if (length(unique(result$samples$reagent.id[result$samples$well.number %in% 
        ladder.wells])) == 1 && "Electronic Ladder" %in% result$samples$sample.name[result$samples$well.number %in% 
        ladder.wells]) {
        result$samples$ladder.well <- ladder.wells[1]
    }
    else if (all.equal(as.character(result$samples$reagent.id[result$samples$well.number %in% 
        ladder.wells]), levels(result$samples$reagent.id))) {
        for (ladder.well in ladder.wells) result$samples$ladder.well[result$samples$reagent.id == 
            result$samples$reagent.id[result$samples$well.number == 
                ladder.well]] <- ladder.well
    }
    result
}

# Function: annotate.electrophoresis
annotate.electrophoresis <- function (electrophoresis, annotations, header = TRUE, row.names = NULL, 
    sep = "\t", stringsAsFactors = TRUE, ...) 
{
    if (any(class(annotations) %in% c("character", "connection"))) 
        annotations <- read.table(annotations, header = header, 
            row.names = row.names, sep = sep, stringsAsFactors = F, 
            ...)
    stopifnot(`empty annotations` = ncol(annotations) > 1, `no label recognized in annotations` = colnames(annotations)[1] %in% 
        colnames(electrophoresis$samples), `duplicate annotations` = anyDuplicated(annotations[, 
        1]) == 0)
    identifiers <- as.character(electrophoresis$samples[, colnames(annotations)[1]])
    for (col in 2:ncol(annotations)) {
        annotation.lookup <- annotations[, col]
        if (stringsAsFactors && class(annotation.lookup) == "character") 
            annotation.lookup <- factor(annotation.lookup, levels = unique(annotation.lookup))
        names(annotation.lookup) <- annotations[, 1]
        electrophoresis$samples[, colnames(annotations)[col]] <- annotation.lookup[identifiers]
    }
    electrophoresis
}

# Function: integrate.custom
integrate.custom <- function (electrophoresis, lower.bound = -Inf, upper.bound = Inf, 
    bound.variable = "length", sum.variable = "molarity") 
as.vector(by(electrophoresis$data, electrophoresis$data$sample.index, 
    function(data.subset) {
        in.this.region <- in.custom.region(data.subset, lower.bound, 
            upper.bound, bound.variable)
        if (sum(in.this.region) == 0) 
            NA
        else sum(data.subset[[sum.variable]][in.this.region])
    }))

# Function: qplot.electrophoresis
qplot.electrophoresis <- function (electrophoresis, x = "length", y = "molarity", ..., 
    log = "", normalize = c("none", "total", "window"), facets = ~sample.index, 
    margins = FALSE, scales = "fixed", geom = c("line", "area"), 
    include.ladder = FALSE, include.markers = FALSE, lower.marker.spread = 10, 
    xlim = c(NA, NA), ylim = c(NA, NA), show.peaks = c("all", 
        "markers", "none"), region.alpha = 0.2, area.alpha = 0.2, 
    title = NULL, xlab = NULL, ylab = NULL) 
{
    normalize <- match.arg(normalize)
    geom <- match.arg(geom)
    show.peaks <- match.arg(show.peaks)
    if (!include.ladder) 
        electrophoresis <- subset(electrophoresis, well.number != 
            ladder.well)
    if (!include.markers) 
        electrophoresis$data <- electrophoresis$data[which(between.markers(electrophoresis, 
            lower.marker.spread)), ]
    electrophoresis$data <- cbind(electrophoresis$data, electrophoresis$samples[electrophoresis$data$sample.index, 
        ])
    if (!is.null(electrophoresis$peaks)) 
        electrophoresis$peaks <- cbind(electrophoresis$peaks, 
            electrophoresis$samples[electrophoresis$peaks$sample.index, 
                ])
    if (!is.null(electrophoresis$regions)) 
        electrophoresis$regions <- cbind(electrophoresis$regions, 
            electrophoresis$samples[electrophoresis$regions$sample.index, 
                ])
    electrophoresis$data$x.value <- electrophoresis$data[[x]]
    if (normalize %in% c("none", "window")) 
        electrophoresis$data <- electrophoresis$data <- electrophoresis$data[which((!is.na(electrophoresis$data$x.value)) & 
            (is.na(xlim[1]) | electrophoresis$data$x.value >= 
                xlim[1]) & (is.na(xlim[2]) | electrophoresis$data$x.value <= 
            xlim[2])), ]
    electrophoresis$data$y.normalized <- if (normalize %in% c("total", 
        "window")) 
        normalize.proportion(electrophoresis, y, lower.marker.spread)
    else electrophoresis$data[[y]]
    electrophoresis$data$y.scaled <- if (y == "fluorescence") 
        electrophoresis$data$y.normalized
    else differential.scale(electrophoresis, x, "y.normalized")
    if (normalize == "total") 
        electrophoresis$data <- electrophoresis$data <- electrophoresis$data[which((!is.na(electrophoresis$data$x.value)) & 
            (is.na(xlim[1]) | electrophoresis$data$x.value >= 
                xlim[1]) & (is.na(xlim[2]) | electrophoresis$data$x.value <= 
            xlim[2])), ]
    electrophoresis$data <- electrophoresis$data[which((!is.na(electrophoresis$data$y.scaled)) & 
        ((!log %in% c("x", "xy")) | electrophoresis$data$x.value > 
            0) & ((!log %in% c("y", "xy")) | electrophoresis$data$y.scaled > 
        0)), ]
    this.plot <- ggplot(electrophoresis$data)
    if (!is.null(facets) & !is.na(region.alpha) & !is.null(electrophoresis$regions)) 
        this.plot <- this.plot + geom_rect(aes_(xmin = as.name(paste0("lower.", 
            x)), xmax = as.name(paste0("upper.", x)), ymin = -Inf, 
            ymax = Inf), data = electrophoresis$regions, alpha = region.alpha)
    this.plot <- this.plot + switch(geom, line = geom_line(if (is.null(facets) && 
        !any(c("color", "colour") %in% names(eval(substitute(alist(...)))))) aes(x = x.value, 
        y = y.scaled, group = sample.index, color = sample.name, 
        ...) else aes(x = x.value, y = y.scaled, group = sample.index, 
        ...)), area = if (is.null(facets) && !"fill" %in% names(eval(substitute(alist(...))))) geom_area(aes(x = x.value, 
        y = y.scaled, group = sample.index, fill = sample.name, 
        ...), alpha = area.alpha) else geom_area(aes(x = x.value, 
        y = y.scaled, group = sample.index, ...), alpha = area.alpha))
    if (!is.null(facets) & show.peaks != "none" & !is.null(electrophoresis$peaks)) {
        peak.data <- subset(cbind(electrophoresis$data, peak = in.peaks(electrophoresis)), 
            !is.na(peak))
        peak.data$peak.observations <- electrophoresis$peaks$peak.observations[peak.data$peak]
        if (show.peaks == "markers") 
            peak.data <- subset(peak.data, peak.observations %in% 
                union(LOWER.MARKER.NAMES, UPPER.MARKER.NAMES))
        this.plot <- this.plot + if (length(unique(peak.data$peak.observations)) > 
            1) 
            geom_area(aes(x = x.value, y = y.scaled, group = peak, 
                fill = peak.observations), data = peak.data)
        else geom_area(aes(x = x.value, y = y.scaled, group = peak), 
            data = peak.data, fill = "darkgray")
    }
    if (!is.null(facets)) 
        this.plot <- this.plot + if (length(facets) == 2) 
            facet_wrap(facets, scales = scales, labeller = if (facets == 
                ~sample.index) 
                labeller.electrophoresis(electrophoresis)
            else "label_value")
        else facet_grid(facets, margins = margins, scales = scales, 
            labeller = if (facets == (. ~ sample.index) || facets == 
                (sample.index ~ .)) 
                labeller.electrophoresis(electrophoresis)
            else "label_value")
    if (!all(is.na(xlim))) 
        this.plot <- this.plot + lims(x = xlim)
    if (!all(is.na(ylim))) 
        this.plot <- this.plot + lims(y = ylim)
    if (log %in% c("x", "xy")) 
        this.plot <- this.plot + scale_x_log10()
    if (log %in% c("y", "xy")) 
        this.plot <- this.plot + scale_y_log10()
    this.plot <- this.plot + labs(x = if (!is.null(xlab)) 
        xlab
    else variable.label(electrophoresis, x), y = if (!is.null(ylab)) 
        ylab
    else variable.label(electrophoresis, (if (normalize != "none") 
        paste("proportion of", normalize, y)
    else y), if (y == "fluorescence") 
        NULL
    else x), title = title)
    if (x %in% c("distance", "relative.distance")) 
        this.plot <- this.plot + scale_x_reverse()
    this.plot
}

# Function: between.markers
between.markers <- function (electrophoresis, lower.marker.spread = 10) 
{
    result <- rep(F, nrow(electrophoresis$data))
    for (lower.marker in which(electrophoresis$peaks$peak.observations %in% 
        LOWER.MARKER.NAMES)) {
        lower.bound <- if (is.na(electrophoresis$peaks$upper.length[lower.marker])) 
            electrophoresis$peaks$length[lower.marker]
        else lower.marker.spread * (electrophoresis$peaks$upper.length[lower.marker] - 
            electrophoresis$peaks$length[lower.marker])
        result[electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[lower.marker] & 
            electrophoresis$data$length > lower.bound] <- T
    }
    for (upper.marker in which(electrophoresis$peaks$peak.observations %in% 
        UPPER.MARKER.NAMES)) {
        upper.bound <- if (is.na(electrophoresis$peaks$lower.length[upper.marker])) 
            electrophoresis$peaks$length[upper.marker]
        else electrophoresis$peaks$lower.length[upper.marker]
        result[electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[upper.marker] & 
            electrophoresis$data$length >= upper.bound] <- F
    }
    result
}

# Function: differential.scale
differential.scale <- function (electrophoresis, x, y) 
{
    stopifnot(`sample indexes out of order` = all(diff(electrophoresis$data$sample.index) >= 
        0))
    delta.x <- unlist(by(electrophoresis$data, electrophoresis$data$sample.index, 
        function(data.subset) c(NA, diff(data.subset[[x]])), 
        simplify = F))
    if (all(delta.x < 0, na.rm = T)) 
        delta.x <- -delta.x
    else stopifnot(`x-values out of order` = all(delta.x > 0, 
        na.rm = T))
    electrophoresis$data[[y]]/delta.x
}

# Function: in.peaks
in.peaks <- function (electrophoresis) 
{
    result <- rep(NA, nrow(electrophoresis$data))
    if (!is.null(electrophoresis$peaks)) 
        for (i in seq(nrow(electrophoresis$peaks))) result[which(in.peak(electrophoresis, 
            i))] <- i
    result
}

# Function: labeller.electrophoresis
labeller.electrophoresis <- function (electrophoresis) 
function(factor.frame) {
    stopifnot(`multiple columns in factor frame` = ncol(factor.frame) == 
        1, `non-integer factor frame` = class(factor.frame[, 
        1]) == "integer")
    list(as.character(electrophoresis$samples$sample.name[factor.frame[, 
        1]]))
}

# Function: normalize.proportion
normalize.proportion <- function (electrophoresis, variable, lower.marker.spread = 5) 
{
    which.usable <- which(between.markers(electrophoresis) & 
        !is.na(electrophoresis$data[[variable]]))
    subset.usable <- electrophoresis$data[which.usable, ]
    sample.sums <- as.vector(by(subset.usable, subset.usable$sample.index, 
        function(data.subset) sum(data.subset[[variable]])))
    result <- rep(NA, nrow(electrophoresis$data))
    result[which.usable] <- subset.usable[[variable]]/sample.sums[subset.usable$sample.index]
    result
}

# Function: variable.label
variable.label <- function (electrophoresis, variable, variable2 = NULL) 
if (is.null(variable2)) switch(variable, time = "time (s)", aligned.time = "aligned time relative to markers (s)", 
    distance = "distance migrated", relative.distance = "distance migrated relative to markers", 
    fluorescence = "fluorescence", length = {
        length.units <- unique(sapply(electrophoresis$assay.info, 
            function(x) x$length.unit))
        if (length(length.units) == 1) paste0("length (", length.units, 
            ")") else {
            warning("incompatible length units")
            "length"
        }
    }, concentration = {
        concentration.units <- unique(sapply(electrophoresis$assay.info, 
            function(x) x$concentration.unit))
        if (length(concentration.units) == 1) paste0("concentration (", 
            concentration.units, ")") else {
            warning("incompatible concentration units")
            "concentration"
        }
    }, molarity = {
        molarity.units <- unique(sapply(electrophoresis$assay.info, 
            function(x) x$molarity.unit))
        if (length(molarity.units) == 1) paste0("molarity (", 
            molarity.units, ")") else {
            warning("incompatible molarity units")
            "molarity"
        }
    }, variable) else paste(variable.label(electrophoresis, variable), 
    "per", variable.label(electrophoresis, variable2))
